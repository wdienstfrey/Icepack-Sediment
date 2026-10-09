"""Check a candidate straight corridor against BedMachine; do not delineate a basin.

The approximate 1953 front is inferred from a modern flowline with a known fjord
extension and the configured advance offset. It is a geometry probe, not a
mapped historical boundary. Inputs remain outside git.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.plotting import style
from src.terminus import chronology

parser=argparse.ArgumentParser()
parser.add_argument("bedmachine")
parser.add_argument("flowline", help="CSV with distance_km,x,y, oriented downstream")
parser.add_argument("--fjord-extension-m",type=float,default=2000)
parser.add_argument("--output",default="results/regional_geometry")
args=parser.parse_args()
line=np.genfromtxt(args.flowline,delimiter=",",names=True)
s=line['distance_km']*1000
if np.any(np.diff(s)<=0): raise ValueError("Flowline must increase downstream")
front_s=s[-1]-args.fjord_extension_m-chronology()[1][-1]
if front_s<30000: raise ValueError("Flowline too short for the 30 km corridor")
center=lambda t: np.array([np.interp(t,s,line['x']),np.interp(t,s,line['y'])])
front,upstream=center(front_s),center(front_s-30000)
direction=(front-upstream)/np.linalg.norm(front-upstream)
normal=np.array([-direction[1],direction[0]])
origin=front-30000*direction-1500*normal
x=np.linspace(0,31820,214);y=np.linspace(0,3000,21)
X,Y=np.meshgrid(x,y,indexing='ij')
points=origin+X[...,None]*direction+Y[...,None]*normal
xmin,ymin=points.reshape(-1,2).min(axis=0)-1000
xmax,ymax=points.reshape(-1,2).max(axis=0)+1000
with xr.open_dataset(args.bedmachine) as dataset:
    yslice=slice(ymax,ymin) if dataset.y.values[0]>dataset.y.values[-1] else slice(ymin,ymax)
    crop=dataset[['mask','thickness','surface','bed']].sel(x=slice(xmin,xmax),y=yslice).load()
xy={'x':xr.DataArray(points[...,0].ravel(),dims='point'),
    'y':xr.DataArray(points[...,1].ravel(),dims='point')}
mask=crop['mask'].interp(xy,method='nearest').values.reshape(X.shape)
h=crop['thickness'].interp(xy).values.reshape(X.shape)
initial=X<=30000
ice=np.isin(mask,[2,3])
summary={'status':'geometry probe only; not prepared model input',
    'bedmachine':str(Path(args.bedmachine)),
    'flowline_sha256':hashlib.sha256(Path(args.flowline).read_bytes()).hexdigest(),
    'assumption':'Modern flowline minus 2 km fjord extension and 1820 m approximate advance',
    'length_1953_m':30000,'width_m':3000,'origin_easting':float(origin[0]),
    'origin_northing':float(origin[1]),
    'downstream_azimuth_deg_from_east':float(np.rad2deg(np.arctan2(direction[1],direction[0]))),
    'initial_nonice_fraction':float(np.mean(~ice[initial])),
    'initial_nonpositive_thickness_fraction':float(np.mean(h[initial]<=0)),
    'initial_missing_fraction':float(np.mean(~np.isfinite(h[initial]))),
    'initial_min_thickness_m':float(np.nanmin(h[initial]))}
out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
style()
fig,ax=plt.subplots(figsize=(6,6),layout='constrained')
ax.pcolormesh(crop.x/1000,crop.y/1000,crop['mask'],vmin=0,vmax=4,cmap='Greys',shading='auto')
ax.plot(line['x']/1000,line['y']/1000,color='#56B4E9',label='Modern flowline')
outline=np.array([origin,origin+30000*direction,origin+30000*direction+3000*normal,
                  origin+3000*normal,origin])
ax.plot(outline[:,0]/1000,outline[:,1]/1000,color='#E69F00',label='Candidate initial corridor')
ax.set(xlim=(xmin/1000,xmax/1000),ylim=(ymin/1000,ymax/1000),aspect='equal',
       xlabel='EPSG:3413 easting (km)',ylabel='EPSG:3413 northing (km)',
       title='Straight-corridor geometry probe')
ax.legend(fontsize=8)
fig.savefig(out/'candidate_corridor.pdf');fig.savefig(out/'candidate_corridor.png')
print(json.dumps(summary,indent=2))
