#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand, sqrt
from fractions import Fraction
import sys

"""
We are analyzing the degree=3 (1,1,1,0,0) Ro3-invariant rule with 7 QPs.
"""

# Create sympy objects for the dofs of this rule.
wc, wv, we, xe = sympify('wc, wv, we, xe')

# Create sympy for the equations. They are in the form "eqn = 0".
eqns = []
# (0,0)
eqns.append(wc + 3*wv + 3*we - Fraction(1,2))
# (2,0)
eqns.append(Fraction(1,9)*wc + wv + we*(2*xe**2 - 2*xe + 1) - Fraction(1,12))
# (3,0)
eqns.append(Fraction(1,27)*wc + wv + we*(3*xe**2 - 3*xe + 1) - Fraction(1,20))
# (2,1)
eqns.append(Fraction(1,27)*wc + we*(xe**3 - 2*xe**2 + xe) - Fraction(1,60))

# Print to verify inputs.
for eqn in eqns:
    print('{}'.format(eqn))

# Verify that a solution we found by hand satisfies the original equations
print('---')
print('Checking solution [wc=9/40, wv=1/40, we=1/15, xe=1/2]')
for eqn in eqns:
    verified = eqn.subs([(wc, Fraction(9,40)),
                         (wv, Fraction(1,40)),
                         (we, Fraction(1,15)),
                         (xe, Fraction(1,2))])
    print('verified = {}, should be 0.'.format(verified))

# Take linear combination of equations to cancel all variables except we, xe
res = simplify(eqns[1] - eqns[2] - 2*eqns[3])

print('---')
print('res={}'.format(res))

# At this point, we have shown that there are two possible solutions,
# xe = (1/2, 1) which corresponds to either QPs at the midpoints of
# each edge, or a "degenerate" edge orbit with points at the vertices.
# We reject the solution xe=1 in this case since we already have a
# a vertex orbit.

# Eliminate wc from eqns (2) and (3)
eqns[1] = simplify(eqns[1] - Fraction(1,9)*eqns[0])
eqns[2] = simplify(eqns[2] - Fraction(1,27)*eqns[0])

print('---')
for eqn in eqns:
    print('{}'.format(eqn))

# Eliminate wv by taking a linear combination of (2) and (3)
res = simplify(Fraction(4,3) * eqns[1] - eqns[2])

print('---')
print('res={}'.format(res))

# Now substitute in xe=1/2 to show that we = 1/15:
res = simplify(res.subs(xe, Fraction(1,2)))

print('---')
print('res={}'.format(res))

# Substitute xe=1/2, we=1/15 into eqn (2) to get wv=1/40
res = simplify(eqns[1].subs([(we, Fraction(1,15)), (xe, Fraction(1,2))]))

print('---')
print('res={}'.format(res))

# Substitute xe=1/2, we=1/15, wv=1/40 into eqn (1) to get wc=9/40
res = simplify(eqns[0].subs([(we, Fraction(1,15)), (wv, Fraction(1,40))]))

print('---')
print('res={}'.format(res))
