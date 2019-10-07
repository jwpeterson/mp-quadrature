#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand, sqrt
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
print('res = {}'.format(res))

# Eliminate w2
#res = simplify(res - w1 - w2 + Fraction(1,6))
res = simplify(res.subs(w2, Fraction(1,6) - w1))
print('res = {}'.format(res))
res = expand(res)
print('res = {}'.format(res))
res = collect(res, w1)
print('res = {}'.format(res))
