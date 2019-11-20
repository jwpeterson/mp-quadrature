#!/usr/bin/env python
import sympy
from sympy.abc import a, x, t
from fractions import Fraction
import sys

# a = alpha
# s = sigma
s = (6*a**2 - 4*a + Fraction(1,2)) / (a * (15*a**2 -9*a  + 1))

# Cubic polynomial
# x = x1
p = 15*s*x**3 - (9*s + 6)*x**2 + (s+4)*x - Fraction(1,2)

# Divide out the factor x-a
result_tuple = sympy.pdiv(p, x-a)
result = sympy.simplify(result_tuple[0])
remainder = sympy.simplify(result_tuple[1])

# Confirm that the remainder is zero
if (remainder != 0):
    print('remainder is nonzero')
    sys.exit(1)

# Print the result, collected in x
print('result={}'.format(sympy.collect(result, x)))

# The denominator is the same as the denominator of sigma, which we assumed
# to be nonzero. Therefore, we can just worry about the numerator when searching
# for the roots.
# x**2*(180*a**2 - 120*a + 15) + x*(-120*a**2 + 75*a - 9) + (15*a**2 - 9*a + 1)
numerator = sympy.numer(result)
print('numerator={}'.format(sympy.collect(numerator, x)))

# roo is a dict, and it seems that the keys are the roots and the dict
# values are their multiplicities (maybe?) Also, note there are two
# roots, so this should correspond to the possibility of having
# multiple solutions for a given value of alpha which we observed
# graphically.
# 1.) -sqrt(3)*sqrt((5*a - 1)*(240*a**3 - 240*a**2 + 75*a - 7))/(30*(12*a**2 - 8*a + 1)) + (40*a**2 - 25*a + 3)/(10*(2*a - 1)*(6*a - 1))
# 2.) sqrt(3)*sqrt((5*a - 1)*(240*a**3 - 240*a**2 + 75*a - 7))/(30*(12*a**2 - 8*a + 1)) + (40*a**2 - 25*a + 3)/(10*(2*a - 1)*(6*a - 1))
root_dict = sympy.simplify(sympy.roots(numerator, x))
for r in root_dict:
    print(r)
    # Substitute in some specific values.
    # x1_onefifth = r.subs([(a,Fraction(1,5))])
    # print('x1(alpha=1/5)={}'.format(x1_onefifth))

    # x1_rootA = sympy.simplify(r.subs([(a, (25 + sympy.sqrt(145))/80)]))
    # print('x1(alpha=root of A)={}, ~ {}'.format(x1_rootA, x1_rootA.evalf()))

# When alpha = alpha3 = (2 + cos(theta/3)) / 6, we would like to know the analytical
# value of x1 which should match what we called "min_alpha" in early versions of the
# analysis. In the formula below, "t" is theta. The discriminant Delta(alpha) = 0
# because alpha1 is one of the roots of the cubic polynomial in Delta.
# Note: theta = pi - arctan(3,4) ~ 2.498091544796508851659834154562

# alpha3=cos(t/3)/6 + 1/3
# x1_alpha3=(-cos(t/3)/12 - cos(2*t/3)/6 + 1/10)/sin(t/3)**2
# x1_alpha3_numerical=0.109039009072877
print('---')
alpha3 = (2 + sympy.cos(t/3)) / 6
print('alpha3={}'.format(alpha3))
x1_alpha3 = sympy.simplify((120*alpha3**2 - 75*alpha3 + 9) / (30 * (2*alpha3 - 1) * (6*alpha3 - 1)))
print('x1_alpha3={}'.format(x1_alpha3))
x1_alpha3_numerical = x1_alpha3.subs([(t,2.498091544796508851659834154562)]).evalf()
print('x1_alpha3_numerical={}'.format(x1_alpha3_numerical))

# When alpha = alpha1 = (4 - ct3 - np.sqrt(3)*st3) / 12
#                     ~ 1.7048618881295388E-01,
# we would like to know the analytical value of x1. From our plots we
# expect it to be larger than 0.6, i.e.  outside the reference
# element.
# alpha1=-sqrt(3)*sin(t/3)/12 - cos(t/3)/12 + 1/3
# x1_alpha1=(sin(t/3 + pi/6)/6 + cos(2*t/3 + pi/3)/3 + 1/5)/(cos(2*t/3 + pi/3) + 1)
# x1_alpha1_numerical=0.659027622374092
print('---')
alpha1 = (4 - sympy.cos(t/3) - sympy.sqrt(3)*sympy.sin(t/3)) / 12
print('alpha1={}'.format(alpha1))
x1_alpha1 = sympy.simplify((120*alpha1**2 - 75*alpha1 + 9) / (30 * (2*alpha1 - 1) * (6*alpha1 - 1)))
print('x1_alpha1={}'.format(x1_alpha1))
x1_alpha1_numerical = x1_alpha1.subs([(t,2.498091544796508851659834154562)]).evalf()
print('x1_alpha1_numerical={}'.format(x1_alpha1_numerical))
