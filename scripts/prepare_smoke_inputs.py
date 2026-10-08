"""Explicitly synthetic geometry for solver tests; never an EKaS reconstruction."""
from pathlib import Path
import json
import numpy as np
x,y = np.linspace(0,32000,129),np.linspace(0,3000,13)
X,Y = np.meshgrid(x,y,indexing="ij")
bed = -500 + 220*(X/32000)**3
surface = 1000-900*X/32000
h = surface-bed
out = Path("data/processed")
out.mkdir(parents=True,exist_ok=True)
np.savez(out/"synthetic.npz",input_kind=np.array("synthetic"),x=x,y=y,bed=bed,thickness=h,smb=np.zeros_like(X),
         ux=np.full_like(X,1000),uy=np.zeros_like(X),friction=np.full_like(X,0.01))
cfg = json.loads(Path("config.json").read_text())
cfg.update(nx=12,ny=4,spinup_years=0.1,dt_years=0.1,end_year=1953.2)
(out/"smoke_config.json").write_text(json.dumps(cfg,indent=2))
print("Synthetic smoke inputs written; not scientific simulation inputs.")
