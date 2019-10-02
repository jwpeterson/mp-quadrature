#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand
from fractions import Fraction
import sys

"""
We are analyzing the degree=3 (0,0,1,1,0) Ro3-invariant rule with 6 QPs.
This rule seems to have no solutions numerically.
"""

# Create sympy objects for the dofs of this rule.
we, xe, wm, xm = sympify('we, xe, wm, xm')

# Create sympy for the equations. They are in the form "eqn = 0".
eqns = []
# (0,0)
eqns.append(3*we + 3*wm - Fraction(1,2))
# (2,0)
eqns.append(we*(2*xe**2 - 2*xe + 1) + wm*(6*xm**2 - 4*xm + 1) - Fraction(1,12))
# (3,0)
eqns.append(we*(3*xe**2 - 3*xe + 1) + wm*(-6*xm**3 + 12*xm**2 - 6*xm + 1) - Fraction(1,20))
# (2,1)
eqns.append(we*(xe**3 - 2*xe**2 + xe) + wm*(3*xm**3 - 3*xm**2 + xm) - Fraction(1,60))

# Print to verify inputs.
for eqn in eqns:
    print('{}'.format(eqn))

# Add/subtract multiples of equations from each other to eliminate wc contribution.
eqns[2] = simplify(3*eqns[1] - 2*eqns[2])

print('---')
for eqn in eqns:
    print('{}'.format(eqn))
