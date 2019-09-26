#!/usr/bin/env python

import numpy as np
from scipy.optimize import fsolve

# This is a hard-coded solver for the degree=3 Ro3-invariant case with
# two edge orbits. There should be no solution to this system of
# equations, and we can show this by deriving the "transformed" system
# of eta variables as done in test.py
def residual(x):
    r = np.zeros(4)
    for i in xrange(2):
        wi = x[2*i]
        xi = x[2*i+1]
        # The primal system of equations
        # (0,0)
        r[0] += 3*wi
        # (2,0)
        r[1] += wi * (2*xi**2 - 2*xi + 1)
        # (3,0)
        r[2] += wi * (3*xi**2 - 3*xi + 1)
        # (2,1)
        r[3] += wi * xi * (xi**2 - 2*xi + 1)

    r = np.subtract(r, [1./2, 1./12, 1./20, 1./60])
    return r

# Computes the Jacobian analytically
def jacobian(x):
    J = np.zeros((4,4))

    # First row is constant
    J[0,0] = 3
    J[0,2] = 3

    for i in xrange(2):
        wi = x[2*i]
        xi = x[2*i+1]

        # Row 1
        J[1, 2*i] = (2*xi**2 - 2*xi + 1) # deriv wrt wi
        J[1, 2*i+1] = wi * (4*xi - 2) # deriv wrt xi

        # Row 2
        J[2, 2*i] = (3*xi**2 - 3*xi + 1) # deriv wrt wi
        J[2, 2*i+1] = wi * (6*xi* - 3) # deriv wrt xi

        # Row 3
        J[3, 2*i] = (xi**3 - 2*xi**2 + xi) # deriv wrt wi
        J[3, 2*i+1] = wi * (3*xi**2 - 4*xi + 1) # deriv wrt xi

    return J

# Solution vector is in the order: [w1, x1, w2, x2]

# Random initial guess
lower_limits = np.zeros(4);
upper_limits = np.array([1./6, 1.0, 1./6, 1.0])
x = np.random.uniform(low=lower_limits, high=upper_limits)

print('initial guess={}'.format(x))

r = residual(x)
print('initial residual={}'.format(r))

# If iflag==1 if solution found, otherwise, check contents of "mesg"
(sol, infodict, iflag, mesg) = fsolve(residual, x, fprime=jacobian, full_output=True)

if iflag != 1:
    print(mesg)
else:
    # We'll never get here, as this system of equations has no
    # solution, not even a solution with negative weights.

    # The solution (if converged) is in the 0'th entry of the returned tuple
    print('Solution = {}'.format(sol))
    print('Number of function evaluations = {}'.format(infodict['nfev']))
    print('Final residual = {}'.format(infodict['fvec']))
    print('Final residual norm = {}'.format(np.linalg.norm(infodict['fvec'])))
