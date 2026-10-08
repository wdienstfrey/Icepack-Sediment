# Inputs and provenance

No regional raw files are committed. The numerical smoke inputs are synthetic,
and cannot be used as evidence of EKaS agreement. The scientific runs need the
following inputs (all in EPSG:3413; sea level zero):

| Input | Paper source / access | Use and limitations |
|---|---|---|
| Bed and surface | [BedMachine Greenland v5](https://nsidc.org/data/idbmg4/versions/5), Earthdata account required | Export `bed` and `surface` to GeoTIFF; paper uses a corrected bed rather than the stock product. |
| Local bed corrections | Rosier 2025 and Millan et al. 2018, cited in the [paper appendix](https://doi.org/10.1017/jog.2026.10172) | The authors' exact kriged bed was not obtained. A user-supplied corrected bed raster is accepted; stock v5 is an explicit sensitivity/substitute. |
| Surface mass balance | [MAR Greenland archive (producer catalog)](https://orbi.uliege.be/handle/2268/344819), 2021 annual mean | Supply meters ice/year. For kg/m²/year divide by 917; for mm water equivalent/year divide by 917. Hold fixed in both experiments. The catalog links MARv3.14.3/ERA5 data at doi:10.5281/zenodo.19691262. This version is a substitute unless verified against the authors' forcing; extract 2021 and record the version. Do not replace it with a 1990–2020 climatology silently. |
| Observed velocity | [MEaSUREs annual Greenland velocity v4](https://nsidc.org/data/nsidc-0725/versions/4), 2020-12-01 to 2021-11-30 | Raster components in meters/year; rotated to model axes in preprocessing. This initializes velocity and imposes trunk inflow; it is not an inversion. |
| Friction | User-prepared k field (MPa/(m/year)^(1/3)) | Ua optimized slipperiness is not available. Supply a calibrated k raster, or omit it to derive a Weertman driving-stress prior from surface slope, thickness and observed speed. This prior is not an inversion. k=Ua C^(-1/3); do not pass Ua C directly. |
| Front chronology | `terminus.csv`, approximately integrated rates from Figure 13 | 0, 960, 1400, 1820 m at 1953/1985/2007/2021. Approximate, not surveyed coordinates. Replace with digitized/author-provided offsets for closer reproduction. |

NSIDC currently lists BedMachine v5 as retired. If the original files cannot
be obtained, [BedMachine v6](https://nsidc.org/data/idbmg4/versions/6) is the
publicly listed successor and an explicit substitute, not the paper's bed.
Record the version and compare near-front depths before interpreting results.

The exact planform front lines and corrected bed are not bundled. This implementation
uses a straight 30 km by 3 km trunk corridor, not the paper's entire ice-divide
bounded basin. Select origin/azimuth from a georeferenced map so x=30 km is the
1953 front and x increases downstream. Verify the corridor stays within ice/fjord
coverage. Width and orientation are material sensitivity parameters.

Keep downloads under `data/raw/`. Export rasters using GDAL or xarray, respecting
nodata and source coordinate systems. `scripts/prepare_inputs.py` samples them
without filling missing data and records their hashes. Use this manifest schema:

```json
{
  "origin_easting": -40000,
  "origin_northing": -3140000,
  "downstream_azimuth_deg_from_east": -110,
  "maximum_length_m": 32000,
  "width_m": 3000,
  "spacing_m": 150,
  "rasters": {
    "bed": {"path": "data/raw/bed.tif", "scale": 1, "source": "dataset/version/DOI"},
    "surface": {"path": "data/raw/surface.tif", "scale": 1, "source": "dataset/version/DOI"},
    "smb": {"path": "data/raw/mar2021.tif", "scale": 0.0010905125, "source": "MAR kg/m2/year; version required"},
    "ux": {"path": "data/raw/vx.tif", "scale": 1, "source": "MEaSUREs v4 2021"},
    "uy": {"path": "data/raw/vy.tif", "scale": 1, "source": "MEaSUREs v4 2021"},
    "friction": {"path": "data/raw/k.tif", "scale": 1, "source": "document calibration here"}
  }
}
```

**The coordinates above illustrate the schema only; they are not a delineated
EKaS domain.** Set them from verified GIS geometry before preparing scientific
inputs. NPZ fields have shape `(len(x), len(y))`, with x,y increasing, bed in m,
thickness in m, smb in m ice/year, ux/uy in m/year and k in MPa/(m/year)^(1/3).
The preparation computes thickness consistently with hydrostatic flotation.
If no friction raster is given, k=rho_i g H |grad(s)| / |u_obs|^(1/3)
is a first-order Weertman prior; local slope noise and lateral/longitudinal stress
make it imperfect. The script rejects zero observed speed rather than inventing
a cutoff. Record any preprocessing of velocity and surface slopes.
The geometry is nominally 2021 truncated at 1953, not a measured 1953 surface.

For exact reproduction, obtain the author's domain, front coordinates, bed,
MAR forcing and inverted A/C fields; these cannot be recovered from figure
captions alone. No Earthdata credentials are stored or requested by scripts.
