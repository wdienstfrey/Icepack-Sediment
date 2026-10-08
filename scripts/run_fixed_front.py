"""Run the fixed front experiment from the repository root."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.model import cli
if __name__ == "__main__":
    cli(False)
