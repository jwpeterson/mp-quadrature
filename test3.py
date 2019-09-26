#!/usr/bin/env python

import numpy as np
from scipy.optimize import fsolve

# This is the "transformed" nonlinear system of equations which is
# generated from the (degree=3, dim=4) Ro3-invariant equations by
# introducing the variable tranformation:
# eta1 := w1*x1 + w2*x2
# eta2 := w1*x1**2 + w2*x2**2
# eta3 := w1*x1**3 + w2*x2**3
# Overall my impression is that the transformed system of equations is
# harder to solve (takes more iterations, requires initial guess for
# eta3 in addition to the parameters, does not converge to final
# residuals as close to zero) but it is useful in the analysis of the
# problem.
def residual(x, *eta):
    r = np.zeros(4)
    eta1, eta2, eta3 = eta # unpack
    # print('eta1={}'.format(eta1))
    # print('eta2={}'.format(eta2))
    # print('eta3={}'.format(eta3))
    for i in xrange(2):
        wi = x[2*i]
        xi = x[2*i+1]
        r[0] += 3*wi
        r[1] += wi * xi
        r[2] += wi * xi**2
        r[3] += wi * xi**3

    r = np.subtract(r, [1./2, eta1, eta2, eta3])
    return r

# Computes the Jacobian analytically. "eta" is ignored.
def jacobian(x, *eta):
    J = np.zeros((4,4))

    # First row is constant
    J[0,0] = 3
    J[0,2] = 3

    for i in xrange(2):
        wi = x[2*i]
        xi = x[2*i+1]

        # Row 1
        J[1, 2*i] = xi # deriv wrt wi
        J[1, 2*i+1] = wi # deriv wrt xi

        # Row 2
        J[2, 2*i] = xi**2 # deriv wrt wi
        J[2, 2*i+1] = 2 * wi * xi # deriv wrt xi

        # Row 3
        J[3, 2*i] = xi**3 # deriv wrt wi
        J[3, 2*i+1] = 3 * wi * xi**2 # deriv wrt xi

    return J

# Solution vector is in the order: [w1, x1, w2, x2]

# Choose initial guess and 0 <= eta3 <= 1/24 randomly.
# x = np.random.uniform(low=np.zeros(4), high=np.array([1./6, 0.5, 1./6, 0.5]))

# Choose eta3, x1, w1 randomly. Note: we choose x1 and w1 toward the small ends of
# their true ranges so that when we compute x2 using the cube root below it is more
# likely to be positive.
eta3 = np.random.uniform(low=0, high=1./24)
x1 = np.random.uniform(low=0, high = .5 * 0.5)
w1 = np.random.uniform(low=0, high = .5 * 1./6)
# Choose w2 such that the weight equation is satisfied already.
w2 = 1./6 - w1
# Choose x2 consistent with other parameters.
x2 = np.cbrt((eta3 - w1 * x1**3) / w2)
# Set initial guess vector to be used below
x = [w1, x1, x2, x2]

# A case which converged, Solution = [ 0.02869706  0.02506193  0.13796961  0.45076143]
# eta3 = 0.012636860413
# x = [ 0.11865216,  0.15516944,  0.00627746,  0.42250296]

# 2.)
# eta3 = 0.0169694672369
# x = [ 0.08299942,  0.46539328,  0.1373625,   0.15714816]

# 3.)
# eta3 = 0.0214734559275
# x = [ 0.02059125,  0.09070912,  0.00624562,  0.35629304]

# 4.)
# eta3=0.00430425249905
# x=[0.14005297683451184, 0.2084716657413544, 0.4849539527287688, 0.4849539527287688]

# Compute eta1, eta2 given eta3:
eta1 = 1./40 + 3*eta3
eta2 = 1./360 + 2*eta3

# Make tuple of values for conveniently passing to residual, jacobian.
eta = (eta1, eta2, eta3)

# Compute and print initial guess, residual, and eta values.
r = residual(x, *eta)
print('initial guess={}'.format(x))
print('eta1={}'.format(eta1))
print('eta2={}'.format(eta2))
print('eta3={}'.format(eta3))
# print('initial residual={}'.format(r))

# If iflag==1 if solution found, otherwise, check contents of "mesg"
(sol, infodict, iflag, mesg) = fsolve(residual, x, args=eta, fprime=jacobian, full_output=True)

if iflag != 1:
    print(mesg)
else:
    # The solution (if converged) is in the 0'th entry of the returned tuple
    print('Solution = {}'.format(sol))
    print('Number of function evaluations = {}'.format(infodict['nfev']))
    print('Final residual = {}'.format(infodict['fvec']))
    print('Final residual norm = {}'.format(np.linalg.norm(infodict['fvec'])))

    # Compute J analytically at the solution
    # Jac = jacobian(sol)
    # print('Final Jacobian (analytical) = \n{}'.format(Jac))

    # Compute determinant of analytical Jacobian
    # det = np.linalg.det(Jac)
    # condJ = np.linalg.cond(Jac)
    # (lam, v) = np.linalg.eig(Jac)
    # print('det(J) = {}'.format(det))
    # print('cond(J) = {:.10e}'.format(condJ))
    # print('eig(J) = {}'.format(lam))

    # Compute the transformed variables according to their definitions.
    eta3 = sol[0]*sol[1]**3 + sol[2]*sol[3]**3
    eta2 = sol[0]*sol[1]**2 + sol[2]*sol[3]**2
    eta1 = sol[0]*sol[1] + sol[2]*sol[3]
    # print('eta1={}'.format(eta1))
    # print('eta2={}'.format(eta2))
    # print('eta3={}'.format(eta3))

    # These differences should be approximately zero if our linear algebra was correct!
    print('eta1 - (1./40 + 3 * eta3)={}'.format(eta1 - (1./40 + 3 * eta3)))
    print('eta2 - (1./360 + 2 * eta3)={}'.format(eta2 - (1./360 + 2 * eta3)))
