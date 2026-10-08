"""First-order ALE splitting with a fresh mesh per step and identical node ordering.

The corridor length follows the prescribed front. A Jacobian correction preserves
ice volume on transfer. Transport uses physical velocity minus mesh velocity.
This is an approximation requiring timestep and mesh convergence, not Ua remeshing.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from .geometry import load_inputs, mesh_and_spaces, scalar, evaluate
from .forcing import accumulation
from .terminus import front


def save_history(output, history):
    """Preserve diagnostics even if a later solver step fails."""
    import csv
    names = sorted(set().union(*(row.keys() for row in history)))
    with (output/"history.csv").open("w") as stream:
        writer = csv.DictWriter(stream, names)
        writer.writeheader()
        writer.writerows(history)


def run(inputs, config, output, advancing):
    import firedrake as fd
    import icepack
    from .friction import potential
    if fd.COMM_WORLD.size != 1:
        raise RuntimeError("Use one MPI rank: nodal transfer is currently serial")
    data = load_inputs(inputs)
    cfg = json.loads(Path(config).read_text())
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    if (out/"failure.json").exists():
        (out/"failure.json").unlink()  # This directory now records the new attempt.
    length0, width = cfg["length_1953_m"], cfg["width_m"]
    dt = cfg["dt_years"]
    if dt <= 0 or cfg["spinup_years"] < 0:
        raise ValueError("Positive timestep and nonnegative spinup required")
    model = icepack.models.IceStream(friction=potential)
    model.quadrature_degree = lambda **fields: int(cfg.get("quadrature_degree", 2))
    history = []
    mesh, Q, V = mesh_and_spaces(length0, width, cfg["nx"], cfg["ny"])
    h = scalar(data, "thickness", Q)
    u = fd.Function(V)
    u.dat.data[:, 0] = scalar(data, "ux", Q).dat.data_ro
    u.dat.data[:, 1] = scalar(data, "uy", Q).dat.data_ro
    initial_hash = hashlib.sha256(h.dat.data_ro.tobytes()+u.dat.data_ro.tobytes()).hexdigest()
    metadata = dict(config=cfg, advancing=advancing, initial_hash=initial_hash,
                    input_sha256=hashlib.sha256(Path(inputs).read_bytes()).hexdigest(),
                    status="running", input_kind=str(data.get("input_kind", "supplied")),
                    icepack_file=icepack.__file__)
    (out/"metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    elapsed, length = 0., length0
    total = cfg["spinup_years"] + cfg["end_year"] - 1953
    base_h = base_u = None
    grid_x = np.linspace(0.01, length0-0.01, 151)
    grid_y = np.linspace(0.01, width-0.01, 31)
    common = np.array(np.meshgrid(grid_x, grid_y, indexing="ij")).reshape(2, -1).T
    gate = np.column_stack([np.full(101, length0-1000), np.linspace(0.01, width-0.01, 101)])
    next_log = 0.
    while True:
        bed, k = scalar(data, "bed", Q), scalar(data, "friction", Q)
        surface = icepack.compute_surface(thickness=h, bed=bed)
        solver = icepack.solvers.FlowSolver(model, dirichlet_ids=[1], side_wall_ids=[3, 4])
        inflow = fd.Function(V)
        inflow.dat.data[:, 0] = scalar(data, "ux", Q).dat.data_ro
        inflow.dat.data[:, 1] = scalar(data, "uy", Q).dat.data_ro
        fd.DirichletBC(V, inflow, 1).apply(u)
        u = solver.diagnostic_solve(velocity=u, thickness=h, surface=surface,
                                   fluidity=fd.Constant(cfg["fluidity"]), friction=k,
                                   mu=fd.Constant(cfg["mu"]), side_friction=fd.Constant(cfg["side_friction"]))
        if not np.all(np.isfinite(u.dat.data_ro)):
            raise RuntimeError("Nonfinite velocity")
        volume = float(fd.assemble(h * fd.dx))
        gate_h, gate_u = np.asarray(evaluate(h, gate)), np.asarray(evaluate(u, gate))
        year = 1953 + max(0, elapsed-cfg["spinup_years"])
        row = dict(elapsed=elapsed, year=year, front_m=length-length0, volume_m3=volume,
                   min_h=float(h.dat.data_ro.min()), max_h=float(h.dat.data_ro.max()),
                   gate_h_m=float(trapezoid(gate_h, gate[:, 1])/width),
                   gate_speed_myr=float(trapezoid(np.linalg.norm(gate_u,axis=1),gate[:,1])/width),
                   gate_flux_Gtyr=float(917e-12*trapezoid(gate_h*gate_u[:,0],gate[:,1])))
        history.append(row)
        save_history(out, history)
        if len(history) == 1 or elapsed >= next_log:
            print(f"elapsed={elapsed:.2f} yr, historical year={year:.2f}, "
                  f"front={length-length0:.2f} m, min H={row['min_h']:.2f} m", flush=True)
            next_log = elapsed + 1
        if base_h is None and elapsed >= cfg["spinup_years"]-1e-8:
            base_h, base_u = np.asarray(evaluate(h, common)), np.asarray(evaluate(u, common))
            np.savez(out/"initial.npz", points=common, thickness=base_h, velocity=base_u)
        if elapsed >= total-1e-8:
            np.savez(out/"final.npz", points=common, thickness=np.asarray(evaluate(h, common)),
                     surface=np.asarray(evaluate(surface, common)), velocity=np.asarray(evaluate(u, common)))
            break
        step = min(dt, total-elapsed)
        if elapsed < cfg["spinup_years"]:
            step = min(step, cfg["spinup_years"]-elapsed)
        next_year = 1953 + max(0, elapsed+step-cfg["spinup_years"])
        next_length = length0 + float(front(next_year, advancing, cfg["smooth_front"]))
        old_h, old_u = h.dat.data_ro.copy(), u.dat.data_ro.copy()
        old_volume = volume
        mesh, Q, V = mesh_and_spaces(next_length, width, cfg["nx"], cfg["ny"])
        h = fd.Function(Q)
        h.dat.data[:] = old_h * length/next_length
        relative = fd.Function(V)
        relative.dat.data[:] = old_u
        relative.dat.data[:, 0] -= ((next_length-length)/step
                                   * mesh.coordinates.dat.data_ro[:,0]/next_length)
        a = accumulation(data, Q)
        boundary_h = scalar(data, "thickness", Q)
        transport = icepack.solvers.FlowSolver(model, prognostic_solver_type="implicit-euler",
                                                 dirichlet_ids=[1], side_wall_ids=[3,4])
        h = transport.prognostic_solve(step, thickness=h, velocity=relative,
                                       accumulation=a, thickness_inflow=boundary_h)
        if not np.all(np.isfinite(h.dat.data_ro)) or h.dat.data_ro.min() <= 0:
            raise RuntimeError(f"Nonpositive/nonfinite thickness at elapsed={elapsed+step:.3f} yr; "
                               f"min H={h.dat.data_ro.min():.6g} m; reduce dt; no clipping applied")
        n = fd.FacetNormal(mesh)
        # At outflow use solved h; at inflow use the prescribed thickness.
        boundary_state = fd.conditional(fd.inner(relative,n)<0, boundary_h,h)
        budget = float(fd.assemble(a*fd.dx-boundary_state*fd.inner(relative,n)*fd.ds))
        residual = float(fd.assemble(h*fd.dx))-old_volume-step*budget
        history[-1]["mass_residual_m3"] = residual
        history[-1]["smb_minus_boundary_flux_m3yr"] = budget
        save_history(out, history)
        u = fd.Function(V)
        u.dat.data[:] = old_u
        length, elapsed = next_length, elapsed+step
    save_history(out, history)
    metadata["status"] = "executed"
    (out/"metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")


def cli(advancing):
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", default="data/processed/ekas.npz")
    parser.add_argument("--config", default="config.json")
    parser.add_argument("--output", default="results/"+("advancing_front" if advancing else "fixed_front"))
    args = parser.parse_args()
    try:
        run(args.inputs,args.config,args.output,advancing)
    except Exception as error:
        out = Path(args.output)
        out.mkdir(parents=True, exist_ok=True)
        metadata_path = out/"metadata.json"
        if metadata_path.exists():
            metadata = json.loads(metadata_path.read_text())
            metadata["status"] = "failed"
            metadata_path.write_text(json.dumps(metadata,indent=2)+"\n")
        (out/"failure.json").write_text(json.dumps({"error":repr(error),
            "inputs":args.inputs,"config":args.config,"advancing":advancing},indent=2)+"\n")
        raise
