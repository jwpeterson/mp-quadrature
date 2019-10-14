#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand, sqrt, solve
from fractions import Fraction
from scipy.optimize import fsolve
import sys
import math

"""
We are analyzing the degree=3 (0,0,0,2,0) Ro3-invariant rule with two
median orbits and 6 QPs. We initially analyzed this rule numerically
in "test.py", but here the goal is to analyze it using sympy.
"""

def charpoly(xvec, *args):
    sigma = args[0]
    x = xvec[0]
    # print ('x={}, type(x)={}'.format(x, type(x)))
    # print ('sigma={}, type(sigma)={}'.format(sigma, type(sigma)))
    r1 = (9 + math.sqrt(21)) / 30
    r2 = (9 - math.sqrt(21)) / 30
    resid = x * (x-r1) * (x-r2) - sigma
    # print ('r1={}, type(r1)={}'.format(r1, type(r1)))
    # print ('r2={}'.format(r2))
    # print ('resid={}, type(resid)={}'.format(resid, type(resid)))
    return resid

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

# Compute the "eta" variables based on this solution (see also
# test3.py for more information).
x1_soln = (9 + sqrt(21)) / 30
x2_soln = (9 - sqrt(21)) / 30

print('w1 ~ {}'.format(w1_soln.evalf()))
print('w2 ~ {}'.format(w2_soln.evalf()))
print('x1 ~ {}'.format(x1_soln.evalf()))
print('x2 ~ {}'.format(x2_soln.evalf()))

eta1 = simplify(w1_soln*x1_soln + w2_soln*x2_soln)
eta2 = simplify(w1_soln*x1_soln**2 + w2_soln*x2_soln**2)
eta3 = simplify(w1_soln*x1_soln**3 + w2_soln*x2_soln**3)

print('---')
print('eta1 = {}'.format(eta1))
print('eta2 = {}'.format(eta2))
print('eta3 = {}'.format(eta3))

# Check that
# eta1 = 1/40 + 3*eta3
# eta1 = 1/360 + 2*eta3
# in agreement with our original analysis of this problem.
print('eta1 - (1./40 + 3 * eta3)={} (should be zero)'.format(eta1 - (Fraction(1,40) + 3 * eta3)))
print('eta2 - (1./360 + 2 * eta3)={} (should be zero)'.format(eta2 - (Fraction(1,360) + 2 * eta3)))

print('---')
# We should be able to derive other solutions by choosing two values
# somewhat arbitrarily:
# 1.) xi = w2/w1
# 2.) x2_new = r2 + dx
# where r2 = x2_soln is the second root above, by first computing
# the right-hand side:
# sigma := -xi * (r2 + dx) * (r2 - r1 + dx) * dx.
# and then finding the root x of
# x * (x-r1) * (x-r2) - sigma = 0
# which is closest to r2
xi = 1.56
dx = .02
x2_new = float(x2_soln.evalf() + dx)
sigma = float((-xi * x2_new * ((x2_soln - x1_soln).evalf() + dx) * dx))

print ('sigma={}'.format(sigma))
print ('dx={}'.format(dx))

(x1_new, infodict, iflag, mesg) = fsolve(charpoly,
                                      float(x2_soln.evalf()),
                                      args=sigma,
                                      full_output=True)

# Once xi is chosen, w1 and w2 are given immediately in terms of xi:
w1_new = (1. / 6) / (xi+1)
w2_new = (xi / 6) / (xi+1)

if iflag != 1:
    print(mesg)
else:
    print('w1_new = {}, x1_new = {},\nw2_new = {}, x2_new = {}'.format(w1_new, x1_new[0],
                                                                       w2_new, x2_new))

    # Check that this new "solution" satisfies the original equations...
    for eqn in eqns:
        verified = eqn.subs([(w1, w1_new),
                             (x1, x1_new[0]),
                             (w2, w2_new),
                             (x2, x2_new)])
        print('verified = {}, should be 0.'.format(verified))
