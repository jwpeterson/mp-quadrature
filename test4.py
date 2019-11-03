#!/usr/bin/env python
from sympy import sympify, simplify, collect, expand
from fractions import Fraction
import sys

"""
We are analyzing the degree=3 (1,1,0,1,0) Ro3-invariant rule with 7 QPs.
Numerically this system of equation seems to have "many" solutions, and
once I corrected an error in my notes, the analysis also shows there are
infinitely many possible solutions.
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

# Print to verify inputs.
# for eqn in eqns:
#     print('{}'.format(eqn))

# Check the numerical solutions below, which are from my C++ code.
all_vals = [
    [1.98421381988942e-1, 2.46954169251223e-2, 7.58307890785635e-2, 4.91026498357039e-1],
    [2.18217389369796e-1, 2.49179057570057e-2, 6.90096311197288e-2, 4.97549244286278e-1],
    [5.40290404129568e-2, 2.35569673182140e-2, 1.25100019210800e-1, 4.60158568786255e-1],
    [2.24616398129737e-1, 2.49952695495363e-2, 6.67992644072178e-2, 4.99858126753552e-1]
]

for vals in all_vals:
    print('---')
    print('Checking solution (wc={},wv={},wm={},xm={})'.format(vals[0], vals[1], vals[2], vals[3]))
    for eqn in eqns:
        verified = eqn.subs([(wc, vals[0]), (wv, vals[1]), (wm, vals[2]), (xm, vals[3])])
        print('verified = {}, should be 0.'.format(verified))

# Choose alpha such that:
# 0 < alpha <= 1/2; alpha != 1/3

# Notes:
# .) Rule is PB if alpha > (9 + sqrt(21))/30 ~ 0.45275, otherwise,
#    rule is NB. See different cases below for more details.
# .) If alpha=1/2, median orbit conservatively degnerates to edge (midpoint) orbit.
# .) If (9 - sqrt(21))/30 ~ .14725 < alpha < 1/5, then vertex weight
#    *and* median weights are negative.

# Case I: alpha < (9 - sqrt(21))/30 ~ .14725 -> wv negative
# Ex: wc=0.114795918367,wv=-0.0416666666667,wm=0.170068027211,xm=0.1

# Case II: .14725 < alpha < 0.2 -> wc _and_ wv negative
# Ex: wc=-0.142485822306,wv=-0.00462962962963,wm=0.218791570398,xm=0.18

# Case III: alpha = 0.2 -> Reduces "classical" 4 QP case with wc<0
# Ex: wc=-0.28125,wv=0.0,wm=0.260416666667,xm=0.2

# Case IV: 0.2 < alpha < (9 + sqrt(21))/30 ~ 0.45275 -> wc negative
# wc blows up to negative infinity as xm -> 1/3.
# Ex: wc=-7.875,wv=0.0138888888889,wm=2.77777777778,xm=0.3

# Case V: 0.45275 < alpha <= 0.5 -> All weights positive
# Ex: wc=0.15805785124,wv=0.0243055555556,wm=0.0896751606979,xm=0.48
alpha = .2
wm_numerical = 1. / 120 / alpha / (3*alpha - 1)**2
wv_numerical = (5. * alpha - 1) / (120 * alpha)
wc_numerical = 1./2 - 3*wm_numerical - 3*wv_numerical

print('---')
print('Checking solution (wc={},wv={},wm={},xm={})'.format
      (wc_numerical, wv_numerical, wm_numerical, alpha))
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
