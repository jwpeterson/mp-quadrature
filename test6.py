#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand, sqrt, solve
from fractions import Fraction
import sys

"""
We are analyzing the degree=3 (0,0,0,2,0) Ro3-invariant rule with two
median orbits and 6 QPs. We initially analyzed this rule numerically
in "test.py", but here the goal is to analyze it using sympy.
"""

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
