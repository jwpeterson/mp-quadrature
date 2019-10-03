#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand, sqrt
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

# Verify that a solution we found by hand satisfies the original equations
print('---')
print('Checking solution [we=1/60, xe=1/2, wm=3/20, xm=1/6]')
for eqn in eqns:
    verified = eqn.subs([(we, Fraction(1,60)),
                         (xe, Fraction(1,2)),
                         (wm, Fraction(3,20)),
                         (xm, Fraction(1,6))])
    print('verified = {}, should be 0.'.format(verified))

print('---')
print('Checking solution for Case 2a:')
for eqn in eqns:
    verified = eqn.subs([(we, Fraction(1,6) - Fraction(1,12) / (Fraction(13,25) + sqrt(21)/75)),
                         (xe, 1),
                         (wm, Fraction(1,12) / (Fraction(13,25) + sqrt(21)/75)),
                         (xm, (9 + sqrt(21)) / 30)]).evalf()
    print('verified = {}, should be 0.'.format(verified))

# Add/subtract multiples of equations from each other to eliminate variables.
eqns[2] = simplify(3*eqns[1] - 2*eqns[2])
eqns[3] = simplify(4*eqns[3] - eqns[2])
eqns[3] = simplify(eqns[3] + eqns[1])

print('---')
for eqn in eqns:
    print('{}'.format(eqn))

# At this point, we have shown that there are two possible solutions,
# xe = (1/2, 1) which corresponds to either QPs at the midpoints of
# each edge, or a "degenerate" edge orbit with points at the vertices.

# First, consider the case xe = 1/2 and make that substitution
eqns_case1 = []
for eqn in eqns:
    eqns_case1.append(eqn.subs(xe, Fraction(1,2)))

print('---')
print('Case 1')
for eqn in eqns_case1:
    print('{}'.format(eqn))

# Substitute we = 1/6 - wm in remaining equations.
# This loop does not modify the original?
# for eqn in eqns_case1:
#     eqn = eqn.subs(we, Fraction(1,6) - wm)
for i in xrange(len(eqns_case1)):
    eqns_case1[i] = simplify(eqns_case1[i].subs(we, Fraction(1,6) - wm))

print('---')
print('Case 1')
for eqn in eqns_case1:
    print('{}'.format(eqn))

# The solution to eqns_case1[1] are xm = (1/6, 1/2). Since xe=1/2 has already
# been chosen, we select xm = 1/6
for i in xrange(len(eqns_case1)):
    eqns_case1[i] = simplify(eqns_case1[i].subs(xm, Fraction(1,6)))

print('---')
print('Case 1')
for eqn in eqns_case1:
    print('{}'.format(eqn))

# Solving the remaining equations gives wm=3/20, we=1/60. This concludes
# the solution for Case 1.

# Second, consider the case xe = 1 and make that substitution
eqns_case2 = []
for eqn in eqns:
    eqns_case2.append(eqn.subs(xe, 1))

# Substitute in we = 1/6 - wm
for i in xrange(len(eqns_case2)):
    eqns_case2[i] = simplify(eqns_case2[i].subs(we, Fraction(1,6) - wm))

# Now combine the equations to cancel the constant term.
eqns_case2[2] = eqns_case2[2] - Fraction(1,5)*eqns_case2[1]

print('---')
print('Case 2')
for eqn in eqns_case2:
    print('{}'.format(eqn))

# The solutions are
# 2a.) xm = (9 + sqrt(21)) / 30
# 2b.) xm = (9 - sqrt(21)) / 30

# First Consider case 2a
case_2a = simplify(eqns_case2[1].subs(xm, (9 + sqrt(21)) / 30))
print('case_2a = {}'.format(case_2a))
