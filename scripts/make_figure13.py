"""Generate figures only from actual saved simulations, or chronology alone."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.plotting import comparison,terminus_panel
parser=argparse.ArgumentParser()
parser.add_argument("--fixed",default="results/fixed_front")
parser.add_argument("--advancing",default="results/advancing_front")
parser.add_argument("--output",default="results/figures")
parser.add_argument("--terminus-only",action="store_true")
parser.add_argument("--smoke",action="store_true")
args=parser.parse_args()
terminus_panel(args.output)
if not args.terminus_only:
    comparison(args.fixed,args.advancing,args.output,args.smoke)
