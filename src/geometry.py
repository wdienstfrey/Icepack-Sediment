"""A straight local corridor, x downstream, y across-flow; not a full drainage basin."""
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def load_inputs(path):
    with np.load(path, allow_pickle=False) as archive:
        data = {name: archive[name] for name in archive.files}
    required = ["x", "y", "bed", "thickness", "smb", "ux", "uy", "friction"]
    for name in required:
        if name not in data or not np.all(np.isfinite(data[name])):
            raise ValueError(f"Missing/nonfinite input: {name}")
    for name in ["x", "y"]:
        if not np.all(np.diff(data[name]) > 0):
            raise ValueError(f"{name} must increase")
    shape = (len(data["x"]), len(data["y"]))
    for name in required[2:]:
        if data[name].shape != shape:
            raise ValueError(f"{name} must have shape {shape}")
    if np.min(data["thickness"]) <= 0 or np.min(data["friction"]) < 0:
        raise ValueError("Positive thickness and nonnegative friction required")
    return data


def sample(data, name, points):
    return RegularGridInterpolator((data["x"], data["y"]), data[name],
                                   bounds_error=True)(points)


def mesh_and_spaces(length, width, nx, ny):
    import firedrake as fd
    mesh = fd.RectangleMesh(nx, ny, length, width)
    return mesh, fd.FunctionSpace(mesh, "CG", 1), fd.VectorFunctionSpace(mesh, "CG", 1)


def scalar(data, name, space):
    import firedrake as fd
    coordinates = space.mesh().coordinates.dat.data_ro
    field = fd.Function(space, name=name)
    field.dat.data[:] = sample(data, name, coordinates)
    return field


def evaluate(field, points):
    """Evaluate serial CG1 fields using their actual triangle connectivity.

    Avoid Firedrake point-location roundoff on kilometer-scale corridor meshes.
    Linear interpolation on these triangles is the CG1 finite element field.
    """
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    space = field.function_space()
    coordinates = space.mesh().coordinates.dat.data_ro
    triangles = space.cell_node_map().values
    grid = Triangulation(coordinates[:, 0], coordinates[:, 1], triangles)
    values = field.dat.data_ro
    components = values[:, None] if values.ndim == 1 else values
    result = []
    for component in components.T:
        sampled = LinearTriInterpolator(grid, component)(points[:, 0], points[:, 1])
        if np.any(np.ma.getmaskarray(sampled)):
            raise ValueError("Diagnostic points outside the mesh")
        result.append(np.asarray(sampled))
    return result[0] if values.ndim == 1 else np.column_stack(result)
