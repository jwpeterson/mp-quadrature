#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand
from fractions import Fraction

"""
This script uses sympy to verify that some algebra in my notes is correct.
We are analyzing the degree=3 (1,1,0,1,0) Ro3-invariant rule with 7 QPs.
Numerically this system of equation seems to have "many" solutions, but
when I try to solve it by hand so far I only come up with a single "unique"
solution...
"""

# Create sympy objects for the dofs of this rule.
wc, wv, wm, xm = sympify('wc, wv, wm, xm')

# Create sympy for the equations. They are in the form "eqn = 0".
eqns = []
eqns.append(wc + 3*wv + 3*wm - Fraction(1,2))
eqns.append(Fraction(1,9)*wc + wv + wm*(6*xm**2 - 4*xm + 1) - Fraction(1,12))
eqns.append(Fraction(1,27)*wc + wv + wm*(-6*xm**3 + 12*xm**2 - 6*xm + 1) - Fraction(1,20))
eqns.append(Fraction(1,27)*wc + wv + wm*(3*xm**3 - 3*xm**2 + xm) - Fraction(1,60))

for eqn in eqns:
    print('{}'.format(eqn))

# Add/subtract multiples of equations from each other to eliminate wc contribution.
eqns[1] = simplify(eqns[1] - Fraction(1,9)*eqns[0])
eqns[2] = simplify(eqns[2] - Fraction(1,27)*eqns[0])
eqns[3] = simplify(eqns[3] - Fraction(1,27)*eqns[0])

# Add/subtract multiples of equations from each other to eliminate wv contribution.
eqns[2] = simplify(eqns[2] - Fraction(4,3)*eqns[1])
eqns[3] = simplify(eqns[3] - Fraction(4,3)*eqns[1])

# Eliminate the constant term from eqns[3]
eqns[3] = simplify(eqns[3] - 7*eqns[2])

print('---')
for eqn in eqns:
    print('{}'.format(eqn))

# Verify that the solution (-27/96, 0, 25/96, 1/5) satisfies the reduced eqns.
# Note: it is also a solution to the original equations, but apparently it is
# not unique.
print('---')
for eqn in eqns:
    verified = eqn.subs([(wc, Fraction(-27,96)), (wv,0),
                         (wm, Fraction(25,96)), (xm, Fraction(1,5))])
    print('verified = {}, should be 0.'.format(verified))
