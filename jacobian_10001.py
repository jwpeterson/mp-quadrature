#!/usr/bin/env python
import sympy
from sympy.abc import a, w, x, y, z
from fractions import Fraction
from sympy.matrices import *
import sys

ninth = Fraction(1,9)
one27 = Fraction(1,27)
J = Matrix(4, 4,
           [1, 3, 0, 0,
            ninth, x**2 + y**2 + z**2, 2 * w * (x - z), 2 * w * (y - z),
            one27, x**3 + y**3 + z**3, w*(3*x**2 - 3*z**2), w*(3*y**2 - 3*z**2),
            one27, x**2 * y + z**2 * x + y**2 * z, w*(2*x*y + z**2 - 2*x*z - y**2), w*(2*y*z + x**2 - 2*x*z - y**2)])
print('Jacobian = {}'.format(J))

# # Check that it simplifies to the known form when x=y, z=1-2x.
# # Results:
# # [[1, 3, 0, 0],
# # [1/9, 6*x**2 - 4*x + 1, 6*w*x - 2*w, 6*w*x - 2*w],
# # [1/27, -6*x**3 + 12*x**2 - 6*x + 1, -9*w*x**2 + 12*w*x - 3*w, -9*w*x**2 + 12*w*x - 3*w],
# # [1/27, 3*x**3 - 3*x**2 + x, 9*w*x**2 - 6*w*x + w, 0]])
# J_xequaly = sympy.expand(J.subs([(z,1-2*x), (y,x)]))
# # FIXME: It seems that we get a different result when computing the derivative
# # before substituting vs. substituting first and then computing the derivative,
# # so we may not be able to make this comparison correctly...
# # TODO: Need to figure out how to collect/simplify individual entries in Matrix...
# # J_xequaly[1,2] =
# # J_xequaly[1,2].collect(x)
# # J_xequaly.collect(J_xequaly, x) # Can't call collect on Matrix?
# print('---')
# print('Jacobian (x=y) = {}'.format(J_xequaly))

# Compute determinant, do basic simplification.
det = J.det()
print('---')
print('det(J) = {}'.format(det))

det_xy = sympy.simplify(det.subs([(z,1-x-y)]))
print('---')
print('With z=1-x-y: {}'.format(det))

det_xy = sympy.simplify(det_xy / (2 * w**2))
print('---')
print('To simplify the problem, divide by 2w^2 (we note that w=0 makes Jacobian singular):\n{}'.format(det_xy))

det_xy = sympy.factor(det_xy)
print('---')
print('sympy.factor:\n{}'.format(det_xy))
# Amazingly, this works really well and makes the result much simpler to understand. Is this even true?
# Result: (3*x**2 + 3*x*y - 3*x + 3*y**2 - 3*y + 1)**3/9

# On the assumption that (3x-1)*(3y-1) appears in the factored result, what happens if
# we subtract that out. The result is:
# 3*x**2 + 3*x*y - 3*x + 3*y**2 - 3*y + 1 = (3*x-1)*(3*y-1) + (x-y)**2
# Note the plus sign, that's not as convenient as it could be...
# tmp = sympy.simplify((3*x**2 + 3*x*y - 3*x + 3*y**2 - 3*y + 1) - (3*x-1)*(3*y-1))
# tmp = sympy.factor(tmp)
# print('---')
# print('tmp:\n{}'.format(tmp))

# Q: Is it possible this is equal to: (3x-1) * (3y-1) * (3z-1)?
# A: No, we end up with -27*x**2*y + 9*x**2 - 27*x*y**2 + 36*x*y - 9*x + 9*y**2 - 9*y + 2
# which does not match the factored result above.
# tmp = (3*x-1) * (3*y-1) * (3*z-1)
# tmp = sympy.expand(tmp.subs([(z,1-x-y)]))
# print('---')
# print('(3*x-1) * (3*y-1) * (3*z-1) = \n{}'.format(tmp))

# Keep all three of x, y, and z
# det_xyz = sympy.simplify(det / (2 * w**2))
# print('---')
# print('det_xyz scaled by (2 w^2):\n{}'.format(det_xyz))

# For some reason, factoring when x, y, and z are all present does not
# result in a simple expression...
# det_xyz = sympy.factor(det_xyz)
# print('---')
# print('det_xyz.factor:\n{}'.format(det_xyz))

# # Substitute in x=y=a
# det_xya = sympy.simplify(det.subs([(x,a), (y,a)]))
# print('---')
# print('With x=y=a: {}'.format(det_xya))

# # Substitute in x=a, y=(1-a)/2
# det_xya2 = sympy.simplify(det.subs([(x,a), (y,(1-a)/2)]))
# print('---')
# print('With x=a, y=(1-a)/2: {}'.format(det_xya2))

# det_x = sympy.collect(det, x)
# print('---')
# print('Collect in x:\n{}'.format(det_x))
#
# det_y = sympy.collect(det, y)
# print('---')
# print('Collect in y:\n{}'.format(det_y))

# Check whether 3*x-1 is a factor. Use polynomial division.
# If it gives you a non-zero remainder, maybe it is not a factor...
# print('---')
# (answer, remainder) = sympy.pdiv(det, 3*x-1)
# print('Polynomial division by 3*x-1')
# print('answer = {}'.format(answer))
# print('remainder = {}'.format(remainder))

# Check whether x-y is a factor of the determinant. Use polynomial division.
# If it gives you a non-zero remainder, maybe it is not a factor...
# print('---')
# (answer, remainder) = sympy.pdiv(det, x-y)
# print('Polynomial division by x-y')
# print('answer = {}'.format(answer))
# print('remainder = {}'.format(remainder))

# Check whether "2*x + y - 1" is a factor of the determinant. Use polynomial division.
# If it gives you a non-zero remainder, maybe it is not a factor...
# print('---')
# (answer, remainder) = sympy.pdiv(det, 2*x + y - 1)
# print('Polynomial division by 2*x + y - 1')
# print('answer = {}'.format(answer))
# print('remainder = {}'.format(remainder))

# # We know if y=x then determinant should be zero! For some reason I don't get this,
# # so there could be something wrong with the input.
# det_xequaly = sympy.simplify(det.subs([(x,y)]))
# print('---')
# print('With x=y: {}'.format(det_xequaly))
