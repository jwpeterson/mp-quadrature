#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand
from fractions import Fraction
import sys

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
# (0,0)
eqns.append(wc + 3*wv + 3*wm - Fraction(1,2))
# (2,0)
eqns.append(Fraction(1,9)*wc + wv + wm*(6*xm**2 - 4*xm + 1) - Fraction(1,12))
# (3,0)
eqns.append(Fraction(1,27)*wc + wv + wm*(-6*xm**3 + 12*xm**2 - 6*xm + 1) - Fraction(1,20))
# (2,1)
eqns.append(Fraction(1,27)*wc + wm*(3*xm**3 - 3*xm**2 + xm) - Fraction(1,60))

for eqn in eqns:
    print('{}'.format(eqn))

# Verify that some solutions found numerically satisfy the governing equations.
# It seems to me that the solution should be unique, but for some reason my code
# finds multiple solutions, and I think this is not correct... And indeed, either
# my governing equations above are wrong or they are wrong in the code, because
# the numerical solutions from the code don't satisfy these governing equations.
# (Apparently they do integrate exactly all polynomials of the required degree,
# so I'm leaning towards the equations above being wrong in some way, but then
# it is interesting that they *do* lead to a different quadrature rule solution.)
print('---')
for eqn in eqns:
    verified = eqn.subs([(wc, 1.9842138198894232087756355234257e-1),
                         (wv, 2.4695416925122307603625640854424e-2),
                         (wm, 7.5830789078563585437186508364720e-2),
                         (xm, 4.9102649835703920944141032123298e-1)])
    print('verified = {}, should be 0.'.format(verified))

print('---')
for eqn in eqns:
    verified = eqn.subs([(wc, 2.1821738936979620876016536190816e-1),
                         (wv, 2.4917905757005788534116782239244e-2),
                         (wm, 6.9009631119728808545828097124702e-2),
                         (xm, 4.9754924428627855803597412887029e-1)])
    print('verified = {}, should be 0.'.format(verified))

print('---')
for eqn in eqns:
    verified = eqn.subs([(wc, 5.4029040412956874973240568941493e-2),
                         (wv, 2.3556967318214072640266582461888e-2),
                         (wm, 1.2510001921080030236865322789095e-1),
                         (xm, 4.6015856878625571724228393577505e-1)])
    print('verified = {}, should be 0.'.format(verified))

# Recipe for generating an arbitrary valid solution:
# Pick any value 1/5 <= alpha < 1/2 for xm
# Note: alpha != 1/3, since that makes the numerator in wm_numerical 0
alpha = .46
# Compute other values directly, based on alpha
wm_numerical = 1. / 180 / alpha / (6*alpha**2 - 4*alpha + 2./3)
wv_numerical = 1./24 - 1. / 120 / alpha
wc_numerical = 1./2 - 3*wm_numerical - 3*wv_numerical

print('---')
print('Checking solution ({},{},{},{})'.format(wc_numerical, wv_numerical, wm_numerical, alpha))
for eqn in eqns:
    verified = eqn.subs([(wc, wc_numerical),
                         (wv, wv_numerical),
                         (wm, wm_numerical),
                         (xm, alpha)])
    print('verified = {}, should be 0.'.format(verified))


# Add/subtract multiples of equations from each other to eliminate wc contribution.
eqns[1] = simplify(eqns[1] - Fraction(1,9)*eqns[0])
eqns[2] = simplify(eqns[2] - Fraction(1,27)*eqns[0])
eqns[3] = simplify(eqns[3] - Fraction(1,27)*eqns[0])

print('---')
for eqn in eqns:
    print('{}'.format(eqn))

# Add/subtract multiples of equations from each other to eliminate wv contribution.
eqns[2] = simplify(eqns[2] - Fraction(4,3)*eqns[1])
eqns[3] = simplify(eqns[3] + Fraction(1,6)*eqns[1])

print('---')
for eqn in eqns:
    print('{}'.format(eqn))

# At this point, we recognize that eqns[2] and eqns[3] are constant
# multiples of each other! So, the next step will zero out eqns[3].
eqns[3] = simplify(eqns[3] + Fraction(1,2)*eqns[2])

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
