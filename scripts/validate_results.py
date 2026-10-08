"""Fail on unequal starts, invalid chronology, mass drift, or missing response.

Qualitative response is a scientific check, never a enforced model constraint.
Run this on complete experiments; smoke runs use --allow-short.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.terminus import front

parser = argparse.ArgumentParser()
parser.add_argument("--fixed",default="results/fixed_front")
parser.add_argument("--advancing",default="results/advancing_front")
parser.add_argument("--allow-short",action="store_true")
args = parser.parse_args()
paths = [Path(args.fixed),Path(args.advancing)]
meta = [json.loads((p/"metadata.json").read_text()) for p in paths]
assert meta[0]["initial_hash"]==meta[1]["initial_hash"],"Initial states differ"
assert meta[0]["input_sha256"]==meta[1]["input_sha256"]
assert meta[0]["config"]==meta[1]["config"],"Configurations differ"
initial = [np.load(p/"initial.npz") for p in paths]
for key in ["points","thickness","velocity"]:
    np.testing.assert_allclose(initial[0][key],initial[1][key],rtol=1e-10,atol=1e-8)
history = [np.genfromtxt(p/"history.csv",delimiter=",",names=True) for p in paths]
f,a = history
np.testing.assert_allclose(f["elapsed"],a["elapsed"])
np.testing.assert_allclose(f["front_m"],0,atol=1e-8)
np.testing.assert_allclose(a["front_m"],front(a["year"],smooth=meta[1]["config"]["smooth_front"]),atol=1e-7)
for h in history:
    assert np.all(np.isfinite(h["volume_m3"])) and np.all(h["min_h"]>0)
    residual = h["mass_residual_m3"][:-1]
    assert np.max(np.abs(residual)/h["volume_m3"][:-1])<1e-4,"Mass residual too large"
if not args.allow_short:
    assert a["year"][-1]>=2021-1e-6,"Historical experiment incomplete"
    assert a["gate_h_m"][-1]>f["gate_h_m"][-1],"Expected thickening not reproduced"
    assert a["gate_speed_myr"][-1]<f["gate_speed_myr"][-1],"Expected slowdown not reproduced"
print("Start, chronology, positive thickness and mass-budget checks passed.")
print("Qualitative historical response: "+("not assessed (short run)" if args.allow_short else "passed"))
