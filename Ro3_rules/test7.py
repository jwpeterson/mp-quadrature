#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand, sqrt
from fractions import Fraction
import sys

"""
We are analyzing the degree=3 (1,0,0,0,1) Ro3-invariant rule with one
centroid orbit and one general orbit, with 4 totals QPs. This rule has
a known solution where the general orbit degenerates to a median orbit
with wc=-27/96, wg=25/96, x=1/5, y=1/5
"""

# Create sympy objects for the dofs of this rule.
wc, wg, x, y = sympify('wc, wg, x, y')

# Create sympy for the equations. They are in the form "eqn = 0".
eqns = []
# eqns[0], (0,0)
eqns.append(wc + 3*wg - Fraction(1,2))
# eqns[1], (2,0)
eqns.append(Fraction(1,9)*wc + wg*(x**2 + (1-x-y)**2 + y**2) - Fraction(1,12))
# eqns[2], (3,0)
eqns.append(Fraction(1,27)*wc + wg*(x**3 + (1-x-y)**3 + y**3) - Fraction(1,20))
# eqns[3], (2,1)
eqns.append(Fraction(1,27)*wc + wg*(x**2 * y + (1-x-y)**2*x + y**2 * (1-x-y)) - Fraction(1,60))

# Print to verify inputs.
for eqn in eqns:
    print('{}'.format(eqn))

# Verify that the known solution satsifes these equations.
print('---')
for eqn in eqns:
    verified = eqn.subs([(wc, Fraction(-27,96)),
                         (wg, Fraction(25,96)),
                         (x, Fraction(1,5)),
                         (y, Fraction(1,5))])
    print('verified = {}, should be 0.'.format(verified))

# Expand equations to get ideas about how to take linear combinations.
for i in xrange(4):
    eqns[i] = collect(expand(eqns[i]), wg)

print('---')
for eqn in eqns:
    print('{}'.format(eqn))

# By retaining z=1-x-y in the governing equations, the usefulness of
# the following linear combination became apparent:
res = collect(eqns[2] + 2*eqns[3] - eqns[1], wg)
print('---')
print('res={}'.format(res))

# From Wolfram-alpha, we know that the polynomial in (x,y) in "res"
# above can be factored in the following manner:
res2 = (x-y)*(2*x + y - 1)*(x + 2*y - 1)

# We can confirm this by subtracting and ensuring the remainder is 0.
# This means that the general orbit degenerates to a median orbit in
# which x=y. The factors above are not three independent results since
# picking any one of the three implies that the other two are also true.
res3 = simplify(res - wg*res2)
print('---')
print('res3={}, should be 0.'.format(res3))
