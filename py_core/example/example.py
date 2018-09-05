
import sys
sys.path.append("/Users/theoglauch/Documents/PhD/Pyjnu-ECP/py_core/")
from pyjnu import PyRun
from constants import phys_const
import time

t0 = time.time()

PYJNU = PyRun(
    Bfield=0.088,
    delta=230.,
    R=0.19e15,
    d=540.,
    z=0.116,
)

t1 = time.time()
print('Time for intialization {} s'.format(t1-t0))
PYJNU.solve_steady()
print('Time for finding steady state solution {} s'.format(time.time()-t1))
