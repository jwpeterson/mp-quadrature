#!/usr/bin/env python
from sympy import pdiv, simplify
from sympy.abc import a, x
from fractions import Fraction

# Test from website: https://docs.sympy.org/latest/modules/polys/reference.html
# The return value is a (result, remainder) pair.
# result = pdiv(x**2 + 1, 2*x - 4)
# print('result={}'.format(result))

# a = alpha
# s = sigma
s = (6*a**2 - 4*a + Fraction(1,2)) / (a * (15*a**2 -9*a  + 1))
print('s={}'.format(s))

# Cubic polynomial
# x = x1
p = 15*s*x**3 - (9*s + 6)*x**2 + (s+4)*x - Fraction(1,2)
print('p={}'.format(p))

# Divide out the factor x-a
result_tuple = pdiv(p, x-a)
result = simplify(result_tuple[0])
remainder = simplify(result_tuple[1])
print('result={}'.format(result))
print('remainder={}'.format(remainder))
