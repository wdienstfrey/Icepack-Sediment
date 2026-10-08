"""Fast checks independent of Firedrake."""
import unittest
import numpy as np
from src.terminus import chronology, front, advance_rate
from src.friction import traction

class PhysicsTests(unittest.TestCase):
    def test_chronology(self):
        years, positions = chronology()
        np.testing.assert_allclose(front(years),positions)
        grid=np.linspace(1928,2021,1000)
        self.assertTrue(np.all(np.diff(front(grid))>=0))
        np.testing.assert_allclose(front(grid,False),0)
        self.assertEqual(float(front(1928)),0)
        self.assertEqual(float(front(2050)),positions[-1])
        self.assertTrue(np.all(advance_rate(np.linspace(1953,2021,1000))>0))

    def test_mixed_law_limits(self):
        k,pressure,mu=0.01,1.,0.5
        speed=np.logspace(-12,15,200)
        stress=traction(speed,pressure,k,mu)
        self.assertTrue(np.all(np.diff(stress)>0))
        self.assertTrue(np.all(stress<mu*pressure))
        np.testing.assert_allclose(stress[0],k*speed[0]**(1/3),rtol=1e-5)
        np.testing.assert_allclose(stress[-1],mu*pressure,rtol=1e-3)
        np.testing.assert_allclose(traction(speed,0,k),0)

if __name__=="__main__":
    unittest.main()
