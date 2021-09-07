#!/usr/bin/env python3

import sympy as sy
from fractions import Fraction

# Check whether the orbit:
# {(a,b), (1-a-b,a), (b, 1-a-b)}
# gives the same quadrature contribution for some different polynomials

a, b = sy.sympify('a, b')
x = [(a,b), (1-a-b,a), (b, 1-a-b)]

x2 = 0
y2 = 0
x3 = 0
x2y = 0
xy2 = 0
y3 = 0
x4 = 0
x5 = 0
x4y = 0
x3y2 = 0
x2y3 = 0
xy4 = 0
for i in range(3):
    # shorthand notation
    xi = x[i][0]
    yi = x[i][1]
    # quadratic terms
    x2 += xi**2
    y2 += yi**2
    # cubic terms - only 2 of 4 are linearly independent
    x3 += xi**3
    x2y += xi**2 * yi
    xy2 += xi * yi**2
    y3 += yi**3
    # quartic
    x4 += xi**4
    # quintic - only 2 of 5 are linearly independent
    x5 += xi**5
    x4y += xi**4 * yi
    x3y2 += xi**3 * yi**2
    x2y3 += xi**2 * yi**3
    xy4 += xi * yi**4

# print(f'x2={sy.expand(x2)}')
# print(f'y2={sy.expand(y2)}')

# cubic
# print(f'x3={sy.expand(x3)}')
# print(f'x2y={sy.expand(x2y)}')
# print(f'xy2={sy.expand(xy2)}')
# print(f'y3={sy.expand(y3)}')

# We can find a non-trivial linear combination of cubic and quadratic
# polynomials that gives us zero. This means that these 4 polynomials
# are _not_ linearly independent!
print('')
print('Results:')
print('1.) Linear combination of cubic polynomials which is not LI:')
print(f'x**3 + x**2y + xy**2 - x**2 = {sy.expand(x3 + x2y + xy2 - x2)}')
print('')

# quartic
# print(f'x4={sy.expand(x4)}')

# quintic
# print(f'x5={sy.expand(x5)}')
# print(f'x4y={sy.expand(x4y)}')
# print(f'x3y2={sy.expand(x3y2)}')
# print(f'x2y3={sy.expand(x2y3)}')
# print(f'xy4={sy.expand(xy4)}')

print('2.) Linear combination of quintic polynomials which is not LI:')
print(f'2 x**5 + 5 (x**4 y + x**3 y**2) - 5/2 x**4 - 5 x**2 y + 1/2 = {sy.expand(2*x5 + 5*(x4y + x3y2) - Fraction(5,2)*x4 - 5*x2y + Fraction(1,2))}')
