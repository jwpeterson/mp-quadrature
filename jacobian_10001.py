#!/usr/bin/env python
import sympy
from sympy.abc import w, x, y, z
from fractions import Fraction
from sympy.matrices import *
import sys

ninth = Fraction(1,9)
one27 = Fraction(1,27)
J = Matrix(4, 4,
           [1, 3, 0, 0,
            ninth, x**2 + y**2 + z**2, w * (2*x - 2*z), w * (2*y - 2*z),
            one27, x**3 + y**3 + z**3, w*(3*x**2 - 3*z**2), w*(3*y**2 - 3*z**2),
            one27, x**2 * y + z**2 * x + y**2 * z, w*(2*x*y + z**2 - 2*x*z - y**2), w*(2*y*z + x**2 - 2*x*z - y**2)])
print('Jacobian = {}'.format(J))

# Compute determinant, do basic simplification.
det = sympy.simplify(J.det())
print('---')
print('det(J) = {}'.format(det))

det = sympy.simplify(det.subs([(z,1-x-y)]))
print('---')
print('With z=1-x-y: {}'.format(det))

det = sympy.simplify(det / (2 * w**2))
print('---')
print('Divide by 2w^2 (w=0 makes Jacobian singular): {}'.format(det))

# Check whether 3*x-1 is a factor. Use polynomial division.
# If it gives you a non-zero remainder, maybe it is not a factor...
# print('---')
# (answer, remainder) = sympy.pdiv(det, 3*x-1)
# print('Polynomial division by 3*x-1')
# print('answer = {}'.format(answer))
# print('remainder = {}'.format(remainder))

# Check whether x-y is a factor of the determinant. Use polynomial division.
# If it gives you a non-zero remainder, maybe it is not a factor...
print('---')
(answer, remainder) = sympy.pdiv(det, x-y)
print('Polynomial division by x-y')
print('answer = {}'.format(answer))
print('remainder = {}'.format(remainder))

# Check whether "2*x + y - 1" is a factor of the determinant. Use polynomial division.
# If it gives you a non-zero remainder, maybe it is not a factor...
print('---')
(answer, remainder) = sympy.pdiv(det, 2*x + y - 1)
print('Polynomial division by 2*x + y - 1')
print('answer = {}'.format(answer))
print('remainder = {}'.format(remainder))

# We know if y=x then determinant should be zero! For some reason I don't get this,
# so there could be something wrong with the input.
det_xequaly = sympy.simplify(det.subs([(x,y)]))
print('---')
print('With x=y: {}'.format(det_xequaly))
