"""Sample co-registered EPSG:3413 GeoTIFFs into a rotated straight corridor.

Export bed/surface from BedMachine and MAR 2021 SMB/MEaSUREs velocity to
GeoTIFF first, preserving masks. The manifest records explicit unit conversion
and the corridor orientation. No automatic bed smoothing or filling is applied.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import rasterio

parser = argparse.ArgumentParser()
parser.add_argument("manifest", help="JSON; see data/README.md")
parser.add_argument("--output",default="data/processed/ekas.npz")
args = parser.parse_args()
m = json.loads(Path(args.manifest).read_text())
theta = np.deg2rad(m["downstream_azimuth_deg_from_east"])
co,si = np.cos(theta),np.sin(theta)
x = np.arange(0,m["maximum_length_m"]+m["spacing_m"]/2,m["spacing_m"])
y = np.arange(0,m["width_m"]+m["spacing_m"]/2,m["spacing_m"])
X,Y = np.meshgrid(x,y,indexing="ij")
points = np.column_stack([(m["origin_easting"]+co*X-si*Y).ravel(),
                          (m["origin_northing"]+si*X+co*Y).ravel()])
fields,provenance = {},{}
for name in ["bed","surface","smb","ux","uy","friction"]:
    if name == "friction" and name not in m["rasters"]:
        continue
    spec = m["rasters"][name]
    path = Path(spec["path"])
    with rasterio.open(path) as raster:
        if raster.crs.to_epsg()!=3413:
            raise ValueError(f"{name}: EPSG:3413 required")
        values = np.ma.concatenate(list(raster.sample(points,masked=True)))
        if np.any(np.ma.getmaskarray(values)) or not np.all(np.isfinite(values)):
            raise ValueError(f"{name}: corridor extends beyond valid raster coverage")
        fields[name] = np.asarray(values).reshape(X.shape)*spec["scale"]
    provenance[name] = dict(source=spec["source"],path=str(path),scale=spec["scale"],
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest())
# Rotate projected velocity to local downstream / across-flow coordinates.
u,v = fields["ux"].copy(),fields["uy"].copy()
fields["ux"],fields["uy"] = co*u+si*v,-si*u+co*v
# Surface must lie above flotation, using Icepack's hydrostatic convention.
base = np.minimum(fields["bed"],0)
grounded = fields["surface"]-fields["bed"]
floating = fields["surface"]/(1-917/1024)
fields["thickness"] = np.where(fields["surface"]+base*(1024-917)/917 >= 0,
                                 grounded,floating)
if np.any(fields["thickness"]<=0):
    raise ValueError("Nonpositive thickness: check surface/bed and corridor")
if "friction" not in fields:
    # Weertman driving-stress prior, as in the paper, without an inversion.
    sx,sy = np.gradient(fields["surface"],x,y,edge_order=2)
    stress = 917*9.81*1e-6*fields["thickness"]*np.hypot(sx,sy)
    speed = np.hypot(fields["ux"],fields["uy"])
    if np.any(speed<=0):
        raise ValueError("Driving-stress prior requires positive observed speed")
    fields["friction"] = stress/speed**(1/3)
    provenance["friction"] = {"source":"Weertman driving-stress prior; not inverted"}
out = Path(args.output)
out.parent.mkdir(parents=True,exist_ok=True)
np.savez(out,x=x,y=y,**fields)
out.with_suffix(".json").write_text(json.dumps(dict(manifest=m,rasters=provenance),indent=2))
