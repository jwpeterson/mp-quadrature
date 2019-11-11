#!/usr/bin/env python
import sympy
from sympy.abc import a, x
from fractions import Fraction
# import mpmath

# Test from website: https://docs.sympy.org/latest/modules/polys/reference.html
# The return value is a (result, remainder) pair.
# result = sympy.pdiv(x**2 + 1, 2*x - 4)
# print('result={}'.format(result))

# a = alpha
# s = sigma
s = (6*a**2 - 4*a + Fraction(1,2)) / (a * (15*a**2 -9*a  + 1))
# print('s={}'.format(s))

# Cubic polynomial
# x = x1
p = 15*s*x**3 - (9*s + 6)*x**2 + (s+4)*x - Fraction(1,2)
# print('p={}'.format(p))

# Divide out the factor x-a
result_tuple = sympy.pdiv(p, x-a)
result = sympy.simplify(result_tuple[0])
remainder = sympy.simplify(result_tuple[1])
# print('result={}'.format(result))
# print('remainder={}'.format(remainder))

# The denominator is the same as the denominator of sigma, which we assumed
# to be nonzero. Therefore, we can just worry about the numerator when searching
# for the roots.
numerator = sympy.collect(sympy.numer(result), x)
# print('numerator={}'.format(numerator))

# roo is a dict, and it seems that the keys are the roots and the dict
# values are their multiplicities (maybe?) Also, note there are two
# roots, so this should correspond to the possibility of having
# multiple solutions for a given value of alpha which we observed
# graphically.
root_dict = sympy.simplify(sympy.roots(numerator, x))
for r in root_dict:
    print(r)

# print(roo[1])
# print('roo={}'.format(roo))
# help(roo)

# Extra coefficients
# coeff_a =
# poly = sympy.Poly(numerator, x)
# coe = poly.coeffs()
# print('coe={}'.format(coe))
# roo = poly.roots()
# print('roo={}'.format(roo))
# roo = sympy.solve(poly)
