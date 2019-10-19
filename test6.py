#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand, sqrt, solve
from fractions import Fraction
from scipy.optimize import fsolve
import sys
import math
import numpy as np

"""
We are analyzing the degree=3 (0,0,0,2,0) Ro3-invariant rule with two
median orbits and 6 QPs. We initially analyzed this rule numerically
in "test.py", but here the goal is to analyze it using sympy.
"""

# The polynomials f and g which we use to derive additional solutions.
def f(x):
    return 6*x**2 - 4*x + 0.5

def g(x):
    return x * (15*x**2 - 9*x + 1)

# Create sympy objects for the dofs of this rule.
w1, x1, w2, x2 = sympify('w1, x1, w2, x2')

# Create sympy for the equations. They are in the form "eqn = 0".
eqns = []
# eqns[0], (0,0)
eqns.append(3*w1 + 3*w2 - Fraction(1,2))
# eqns[1], (2,0)
eqns.append(w1*(6*x1**2 - 4*x1 + 1) + w2*(6*x2**2 - 4*x2 + 1) - Fraction(1,12))
# eqns[2], (3,0)
eqns.append(w1*(-6*x1**3 + 12*x1**2 - 6*x1 + 1) + w2*(-6*x2**3 + 12*x2**2 - 6*x2 + 1) - Fraction(1,20))
# eqns[3], (2,1)
eqns.append(w1*(3*x1**3 - 3*x1**2 + x1) + w2*(3*x2**3 - 3*x2**2 + x2) - Fraction(1,60))

# Print to verify inputs.
for eqn in eqns:
    print('{}'.format(eqn))

# Confirm that eqns[1] = 2*eqns[3] + eqns[2], i.e. there are only three
# linearly independent equations.
should_be_zero = simplify(eqns[1] - (2*eqns[3] + eqns[2]))
print('should_be_zero = {}'.format(should_be_zero))

# Take linear combination of the two cubic equations, try to cancel
# the constant term.
res = simplify(7*eqns[3] + eqns[2])
print('---')
print('res = {}'.format(res))

# In res, we should now have:
# sum_i w_i x_i (15*x_i**2 - 9*x_i + 1) = 0
# This polynomial has roots:
# (9 + sqrt(21)) / 30
# (9 - sqrt(21)) / 30
# If we pick x1 as the positive root and x2 as the negative root, the
# equation should be satisifed! Note: We also have to eliminate the weights
should_be_zero = simplify(res.subs([(x1, (9 + sqrt(21)) / 30),
                                    (x2, (9 - sqrt(21)) / 30),
                                    (w2, Fraction(1,6) - w1)]))
print('should_be_zero = {}'.format(should_be_zero))

# Now, let's substitute all this into one of the original eqns and simplify
res = simplify(eqns[3].subs([(x1, (9 + sqrt(21)) / 30),
                             (x2, (9 - sqrt(21)) / 30),
                             (w2, Fraction(1,6) - w1)]))

# Finally, we can use sympy to solve the equation res = 0, for w1
w1_soln = solve(res, w1)[0]
w2_soln = Fraction(1,6) - w1_soln

print('---')
# print('res = {}'.format(res))
print('w1 = {}'.format(w1_soln))
print('w2 = {}'.format(w2_soln))
w_ratio = simplify(w2_soln/w1_soln)
# print('Baseline solution, ratio of weights w2/w1 = {} ~ {}'.format(w_ratio, w_ratio.evalf()))

# Compute the "eta" variables based on this solution (see also
# test3.py for more information).
x1_soln = (9 + sqrt(21)) / 30
x2_soln = (9 - sqrt(21)) / 30

print('---')
print('w1 ~ {}'.format(w1_soln.evalf()))
print('w2 ~ {}'.format(w2_soln.evalf()))
print('x1 ~ {}'.format(x1_soln.evalf()))
print('x2 ~ {}'.format(x2_soln.evalf()))

eta1 = simplify(w1_soln*x1_soln + w2_soln*x2_soln)
eta2 = simplify(w1_soln*x1_soln**2 + w2_soln*x2_soln**2)
eta3 = simplify(w1_soln*x1_soln**3 + w2_soln*x2_soln**3)

# print('---')
# print('eta1 = {}'.format(eta1))
# print('eta2 = {}'.format(eta2))
# print('eta3 = {}'.format(eta3))

# Check that
# eta1 = 1/40 + 3*eta3
# eta1 = 1/360 + 2*eta3
# in agreement with our original analysis of this problem.
print('---')
print('eta1 - (1./40 + 3 * eta3)={} (should be zero)'.format(eta1 - (Fraction(1,40) + 3 * eta3)))
print('eta2 - (1./360 + 2 * eta3)={} (should be zero)'.format(eta2 - (Fraction(1,360) + 2 * eta3)))

# Idea to find more solutions.
#
# Rewrite any two equations as:
# a.) w1/w2 = -f(x_2) / f(x_1)
# b.) w1/w2 = -g(x_2) / g(x_1)
# Where f, g are polynomials in their argument, and f(x_1) != 0, g(x_1) != 0
# This seems to always be possible.
#
# Then we can equate a.) and b.) to eliminate the weights and obtain:
# f(x_2)/f(x_1) = g(x_2)/g(x_1)
# Since the system is underdetermined, set x_2 = alpha, where alpha
# is "close" to the x_2 solution we found analytically. Then we have:
# f(alpha)/g(alpha) = f(x_1)/g(x_1)
#
# Let sigma := f(alpha)/g(alpha). Then we have to solve the equation
# f(x_1) - sigma * g(x_1) = 0  (*)
# for x_1. In general (*) will be a cubic polynomial in x_1, so it
# may up to three real roots. We will choose the root which is closest
# to the x_1 value we found analytically.
#
# Then we can use either a.) or b.) above to solve for the ratio
# w1/w2, and finally we can use the fact that w1 + w2 = 1/6 to solve
# for w1 and w2 separately.

# Choose x2=alpha, where x2 is close to the analytical root, x2_soln
pert = np.random.uniform(low=-.022, high=.022)
print('---')
print('pert={}'.format(pert))
alpha = float(x2_soln) + pert

# Compute sigma based on alpha
sigma = f(alpha) / g(alpha)

# Instead of calling a full nonlinear solver, get all the roots of the cubic
# equation, and choose exactly the one that you want. The input to the roots
# function is the polynomial coefficients, starting with the highest power of x.
roots = np.roots([15*sigma, -(9*sigma + 6), (sigma + 4), -0.5])
print('roots={}'.format(roots))

# Find root which is closes to the analytical solution's x1 value
x1_new = min(roots, key=lambda x : abs(x - float(x1_soln)))

# Exit if we failed to find a real root
if not np.isreal(x1_new):
    print('Imaginary root found, try different initial guess.')
    sys.exit(1)

print('Computing arbitrary solution numerically')
x2_new = alpha

# Check for x1 ~ x2
if (np.abs(x1_new - x2_new) < 1.e-9):
    print('Error: nonlinear solver converged to invalid solution, x1=x2')
    sys.exit(1)

# Now solve for the ratio w1/w2
w1_over_w2 = -f(x2_new) / f(x1_new)
w2_over_w1 = 1. / w1_over_w2
# print('w1_over_w2={}'.format(w1_over_w2))
# print('w2_over_w1={}'.format(w2_over_w1))

# Now solve for w1 and w2 separately
w1_new = 1. / 6 / (1 + w2_over_w1)
w2_new = 1. / 6 / (1 + w1_over_w2)

# Print the solution
print('w1 = {}'.format(w1_new))
print('w2 = {}'.format(w2_new))
print('x1 = {}'.format(x1_new))
print('x2 = {}'.format(x2_new))

# Check that this new "solution" satisfies the original equations...
print('---')
for eqn in eqns:
    verified = eqn.subs({w1:w1_new, w2:w2_new, x1:x1_new, x2:x2_new})
    print('verified = {}, should be 0.'.format(verified))
