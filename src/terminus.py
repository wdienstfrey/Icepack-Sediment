"""Editable chronology; approximate rates read from Figure 13, not digitized fronts."""
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator

CHRONOLOGY = Path(__file__).resolve().parents[1] / "data/terminus.csv"


def chronology(path=CHRONOLOGY):
    rows = np.genfromtxt(path, delimiter=",", names=True)
    years, advance = rows["year"], rows["advance_m"]
    if (len(years) < 2 or not np.all(np.diff(years) > 0)
            or not np.all(np.diff(advance) >= 0) or advance[0] != 0):
        raise ValueError("Chronology must start at zero and increase monotonically")
    return years, advance


def front(year, advancing=True, smooth=True, path=CHRONOLOGY):
    years, positions = chronology(path)
    t = np.clip(year, years[0], years[-1])
    if not advancing:
        return np.zeros_like(t, dtype=float)
    if smooth:
        return PchipInterpolator(years, positions)(t)
    return np.interp(t, years, positions)


def advance_rate(year, smooth=True, path=CHRONOLOGY):
    years, positions = chronology(path)
    if smooth:
        return PchipInterpolator(years, positions).derivative()(year)
    index = np.clip(np.searchsorted(years, year, side="right") - 1, 0, len(years)-2)
    return (np.diff(positions) / np.diff(years))[index]
