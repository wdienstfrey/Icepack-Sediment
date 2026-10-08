"""Ua A.2 mixed law expressed as an Icepack dissipation potential.

k = Ua C**(-1/3), in MPa / (m/year)**(1/3), not Ua slipperiness C.
The derivative of the potential gives |tau| = mu*N*k*v**(1/3)/(mu*N+k*v**(1/3)).
16-point Gauss quadrature after a cubic substitution avoids an endpoint singularity.
"""
import numpy as np


def traction(speed, pressure, k, mu=0.5):
    limit = mu * np.maximum(pressure, 0)
    weertman = k * np.asarray(speed)**(1/3)
    return limit * weertman / (limit + weertman + 1e-30)


def potential(**kwargs):
    import firedrake as fd
    from icepack.constants import ice_density, water_density, gravity
    u, h, surface = (kwargs[name] for name in ["velocity", "thickness", "surface"])
    k = kwargs["friction"]
    base = surface - h
    pressure = fd.max_value(ice_density * gravity * h
                           - water_density * gravity * fd.max_value(-base, 0), 0)
    limit = kwargs.get("mu", 0.5) * pressure
    speed = fd.sqrt(fd.inner(u, u) + 1e-6)  # 0.001 m/year regularization
    nodes, weights = np.polynomial.legendre.leggauss(16)
    action = 0
    for z, weight in zip((nodes + 1)/2, weights/2):
        weertman = k * speed**(1/3) * z
        stress = limit * weertman / (limit + weertman + 1e-12)
        action += float(weight) * stress * 3 * speed * float(z)**2
    return action
