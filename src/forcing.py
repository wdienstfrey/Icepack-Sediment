"""Hold the prepared MAR 2021 annual SMB field fixed, in meters ice per year."""
from .geometry import scalar


def accumulation(data, space):
    return scalar(data, "smb", space)
