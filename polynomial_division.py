#!/usr/bin/env python
import sympy
from sympy.abc import a, x
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
    x1_onefifth = r.subs([(a,Fraction(1,5))])
    print('x1(alpha=1/5)={}'.format(x1_onefifth))
    x1_rootA = sympy.simplify(r.subs([(a, (25 + sympy.sqrt(145))/80)]))
    print('x1(alpha=root of A)={}, ~ {}'.format(x1_rootA, x1_rootA.evalf()))
