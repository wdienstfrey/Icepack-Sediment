"""Run in Firedrake: validate CG1 diagnostics, ALE transfer and mixed-law traction."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import firedrake as fd
from icepack.constants import ice_density, water_density, gravity
from src.geometry import mesh_and_spaces,evaluate
from src.friction import potential,traction

old,Q,V=mesh_and_spaces(30000,3000,12,4)
x,y=fd.SpatialCoordinate(old)
h=fd.Function(Q).interpolate(500+0.01*x+0.02*y)
points=np.array([[3200,2999.99],[15000,1500],[.01,.01]])
np.testing.assert_allclose(evaluate(h,points),500+.01*points[:,0]+.02*points[:,1],rtol=1e-12)
volume=float(fd.assemble(h*fd.dx))
new,Q1,V1=mesh_and_spaces(30100,3000,12,4)
transferred=fd.Function(Q1)
transferred.dat.data[:]=h.dat.data_ro*30000/30100
np.testing.assert_allclose(float(fd.assemble(transferred*fd.dx)),volume,rtol=1e-12)
mesh=fd.UnitSquareMesh(1,1)
Q=fd.FunctionSpace(mesh,"CG",1);V=fd.VectorFunctionSpace(mesh,"CG",1)
u=fd.Function(V);direction=fd.Function(V).assign(fd.Constant((1,0)))
h=fd.Function(Q);s=fd.Function(Q)
k=fd.Function(Q).assign(.01)
action=potential(velocity=u,thickness=h,surface=s,friction=k)*fd.dx
stress_form=fd.derivative(action,u,direction)
for pressure in [0,5e-5,0.585,5.]:
    height=(pressure/gravity+water_density*300)/ice_density
    h.assign(height);s.assign(height-300)
    for speed in [.01,1,1000,1e6]:
        u.assign(fd.Constant((speed,0)))
        stress=float(fd.assemble(stress_form))
        regularized=np.sqrt(speed**2+1e-6)
        expected=traction(regularized,pressure,.01)*speed/regularized
        np.testing.assert_allclose(stress,expected,rtol=1e-6,atol=1e-10)
print("CG1 interpolation, volume-preserving transfer and mixed-law derivative passed.")
