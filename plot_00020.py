import numpy as np
import matplotlib.pyplot as plt
import sys
import cmath # sqrt(negative) returns imaginary number instead of nan
from scipy.optimize import minimize
from scipy.optimize import fsolve

from matplotlib import rcParams
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 12
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

"""
For configuration (0,0,0,2,0), plots x1(alpha) vs. alpha
"""

"""
It's convenient to have Python functions for these polynomials.
"""
def f(x):
    return 6*x**2 - 4*x + 0.5

def g(x):
    return x * (15*x**2 - 9*x + 1)

"""
Given alpha, computes x1(alpha) analytically. Returns a pair of
(possibly complex) values, throws an error for input values of
alpha=1/2 or 1/6 since computing x1 by this method is not well-defined
for those values.
"""
def compute_x1_analytical(alpha):
    # Don't accept input values outside the range
    if (alpha <= 0.) or (alpha >= 0.5) or (np.abs(alpha - 1./6) < 1.e-8):
        print('Invalid input, alpha={}; 0 < alpha < 0.5 is required.'.format(alpha))
        raise RuntimeError('Invalid input, alpha=1/6.'.format(alpha))

    # Using the polynomial_division.py code, we found the following roots analytically.
    # Note 1: we intentionally use cmath.sqrt so that we can obtain an imaginary root
    # instead of a NaN!
    # Note 2: If the roots are real-valued, the second one will always be smaller since
    # in that case B > 0
    # Note 3: The notation here matches the discussion in the paper.
    A = 120*alpha**2 - 75*alpha + 9
    B = cmath.sqrt(3 * (5*alpha-1) * (240*alpha**3 - 240*alpha**2 + 75*alpha - 7))
    den = 30 * (2*alpha - 1) * (6*alpha - 1)
    # print('A={}, B={}, den={}'.format(A,B,den))
    analytical_roots = [(A + B)/den, (A - B)/den]
    # print('analytical_roots_v2={}'.format(analytical_roots))
    return analytical_roots

"""
Given x1 and x2, computes the weights w1 and w2
"""
def compute_weights(x1, x2):
    w1_over_w2 = -f(x2) / f(x1)
    # w1_over_w2 = -g(x2) / g(x1) # Same
    w2_over_w1 = 1. / w1_over_w2
    w1 = 1. / 6 / (1 + w2_over_w1)
    w2 = 1. / 6 / (1 + w1_over_w2)
    return w1, w2

"""
Function used with fsolve to find alpha for which w1(alpha) = w2(alpha) = 1/12
"""
def residual(alpha):
    x1 = compute_x1_analytical(alpha)[0].real
    w1, w2 = compute_weights(x1, alpha)
    return w1 - 1./12

"""
Function used with minimize() to find a minimum value of x1. Returns the
first root, "(A+B)/den", throwing an error if the value has a nonzero imaginary component.
"""
def fmin(alpha):
    roots = compute_x1_analytical(alpha)
    if (np.abs(roots[0].imag) > 1.e-6):
        raise RuntimeError('No valid root found for alpha={}, roots={}.'.format(alpha, roots))
    # Real-valued return
    return roots[0].real

################################################################################

# Some useful constants
r1 = (9 + np.sqrt(21))/30 # .45275
r2 = (9 - np.sqrt(21))/30 # .14725

################################################################################

# Testing: are there solutions for 1/6 < alpha < r1? It appears that
# there are mostly one real root (alpha) and two imaginary roots in this interval,
# but we can go a bit outside the interval? I suspect the minimum is the
# min_x1 which we computed previously...
# ---
# Roots of "A"
# alpha = (25 + np.sqrt(145)) / 80 ~ 0.46302 # A = 0, x1= +/- 0.15701394368037236
# alpha = (25 - np.sqrt(145)) / 80 ~ 0.16198 # A = 0, x1= +/- 0.47470687954772173
# ---
# Roots of "B"
# alpha = 0.2 # B=0, x1=1/3 (points at centroid)
# alpha = 0.170486188812954 # "alpha1", B = 0, x1 = 0.65902762237408341 (outside the region)
# alpha = 0.384033315723485 # "alpha2",B = 0, x1 = 0.23193336855303029 (inside the region)
# alpha = 0.445480495463561 # "alpha3", B = 0, x1 = 0.10903900907287808 (inside the region*)
# * The values 0.445480495463561 and 0.10903900907287808 correspond to the (min_alpha, min_x1)
#   values we had found earlier. The analytical expression for these alphas is very complicated
#   but perhaps it is worth writing down since the value seems to be important?

# Roots of "B", trigonometric forms:
# alpha
mu = np.arctan(3./4)
theta = np.pi - mu
tp = theta + np.pi
ct3 = np.cos(theta/3)
st3 = np.sin(theta/3)
cm3 = np.cos(mu/3)
sm3 = np.sin(mu/3)
ctp3 = np.cos(tp/3)
stp3 = np.sin(tp/3)
alpha1_trig = (4 - ct3 - np.sqrt(3)*st3) / 12
alpha2_trig = (4 - ct3 + np.sqrt(3)*st3) / 12
alpha3_trig = (2 + ct3) / 6
x1_alpha1 = (6 + 5*cm3 - 10*np.cos(2*mu/3))/(60 * sm3**2)
x1_alpha2 = (6 + 5*ctp3 - 10*np.cos(2*tp/3))/(60 * stp3**2)
x1_alpha3 = (6 - 5*ct3 - 10*np.cos(2*theta/3))/(60 * st3**2)
print('alpha1_trig={:.16E}'.format(alpha1_trig)) # 1.7048618881295388E-01
print('alpha2_trig={:.16E}'.format(alpha2_trig)) # 3.8403331572348476E-01
print('alpha3_trig={:.16E}'.format(alpha3_trig)) # 4.4548049546356139E-01
print('x1(alpha1)={:.16E}'.format(x1_alpha1)) # 6.5902762237409251E-01
print('x1(alpha2)={:.16E}'.format(x1_alpha2)) # 2.3193336855303059E-01
print('x1(alpha3)={:.16E}'.format(x1_alpha3)) # 1.0903900907287724E-01

# Testing intervals
# .) 0 < alpha < 1/6 = the PI/NI region we already investigated
# alpha = 0.1 # x1 = 0.44562222747978653, 0.11687777252021356
# .) 1/6 < alpha < alpha1
# alpha = 0.169 # x1 = (0.53366013542957336, 1.2479540208069468) (points outside region)
# .) alpha1 < alpha < 0.2
# alpha = 0.18 # x1 = Imaginary, 0.3984375-0.12540624091494318j (IMAG)
# .) 0.2 < alpha < 1/3
# alpha = 0.25 # x1 = (0.21835034190722738, 0.38164965809277263) (NI)
# .) 1/3 < alpha < alpha2
# alpha = 0.375 # x1 = (0.20944949536696106, 0.27055050463303892) (NI)
# .) alpha2 < alpha < alpha3
# alpha = 0.4 # x1 = Imaginary, 0.2142857-0.0412393049j (IMAG)
# .) alpha3 < alpha < 0.5 = the PI region we already investigated
alpha = 0.48 # x1 = -0.737867, 0.1633997                (PI)
# ---
# Other values
# alpha = r2 - .01
# alpha = r2 # error in roots: array must not contain Infs or NaNs
# alpha = 0.445480495467 # min_x1
# alpha = 0.444 imaginary
### try:
###     analytical_roots = compute_x1_analytical(alpha)
###     w1, w2 = compute_weights(analytical_roots[0].real, alpha)
###     print('---')
###     print('w1={}'.format(w1))
###     print('x1={}'.format(x1))
###     print('w2={}'.format(w2))
###     print('x2={}'.format(alpha))
### except Exception as e:
###     print('Exception caught: ' + str(e))
###     # sys.exit(0)

# Early
# sys.exit(0)

################################################################################

# Find the minimum value of x1(alpha). Interestingly, it does not
# seem to coincide with either of the limits alpha=0 or alpha=1/6.
# method='L-BFGS-B'
# method='Newton-CG' # requires Jacobian
# method='CG' # cannot handle constraints or bounds
result = minimize(fmin, 0.11, method='CG', \
                  options={'disp': False,
                           'gtol' : 1.e-10})

# Extract the min x1 value and the alpha where it occurs.
# Note: the value of x1 at the minimum is "alpha3_trig" which is computed exactly above.
min_alpha = result.x[0]
min_x1 = result.fun
# Compute the corresponding weights at the (alpha, x1(alpha)) solution
w1, w2 = compute_weights(min_x1, min_alpha)

# The solution is:
# w1=0.0995222256262
# x1=0.445480495467
# w2=0.0671444410405
# x2=0.109037640529
print('---')
print('Min x1 solution:')
print('w1={}'.format(w1))
print('x1={}'.format(min_x1))
print('w2={}'.format(w2))
print('x2={}'.format(min_alpha))

################################################################################

# Solve for the value of alpha which gives w1 = w2 = 1/12. Initial guess
# is alpha=1/8 based on the graph of w1(alpha). The solution obtained
# numerically is (note: x2=alpha):
# w1=0.0833333333333
# x1=0.446333855871
# w2=0.0833333333333
# x2=0.126484505776
(sol, infodict, iflag, mesg) = fsolve(residual, 1./8, full_output=True)
alpha_weights_equal = 0
if iflag == 1:
    alpha_weights_equal = sol[0]
    # print('Solution = {:.12E}'.format(alpha_weights_equal))
    # print('Number of function evaluations = {}'.format(infodict['nfev']))
    # print('Final residual norm = {}'.format(np.linalg.norm(infodict['fvec'])))

    # Compute and print full solution for this value of alpha
    x1_weights_equal = compute_x1_analytical(alpha_weights_equal)[0].real
    w1, w2 = compute_weights(x1_weights_equal, alpha_weights_equal)
    # Print current result
    print('---')
    print('Equal weights solution:')
    print('w1={}'.format(w1))
    print('x1={}'.format(x1_weights_equal))
    print('w2={}'.format(w2))
    print('x2={}'.format(alpha_weights_equal))
else:
    print('fsolve not converged.')


################################################################################

# This curve is for 0 <= alpha < 1/6, with clustering near 1/6 to
# capture the rapid change near that value.
# We always take the larger root for this plot which is in
# entry [0] of the returned array. Note: we take the real part
# without verifying that the imaginary part is zero or nearly
# zero, but for consistency we should probably do that.
alphas = np.concatenate((np.linspace(1.e-6, 0.15), np.linspace(0.15, 1./6 - 1.e-6)))
w1 = np.zeros(len(alphas))
w2 = np.zeros(len(alphas))
x1 = np.zeros(len(alphas))

for i in xrange(len(alphas)):
    alpha = alphas[i]
    x1[i] = compute_x1_analytical(alpha)[0].real
    w1[i], w2[i] = compute_weights(x1[i], alpha)

################################################################################

# Starting from alpha=0+eps and taking min roots, we can recover a branch
# of NI rules instead. There is a singularity at alpha=min_alpha where
# one set of weights blows up to +infty while the other set blows up to -infty.
# We always take the real part of the smaller root, but we should verify that
# the imaginary part is zero.
min_root_alphas = np.linspace(1.e-6, min_alpha - 1.e-3)
min_root_w1 = np.zeros(len(min_root_alphas))
min_root_w2 = np.zeros(len(min_root_alphas))
min_root_x1 = np.zeros(len(min_root_alphas))

for i in xrange(len(min_root_alphas)):
    alpha = min_root_alphas[i]
    min_root_x1[i] = compute_x1_analytical(alpha)[1].real
    min_root_w1[i], min_root_w2[i] = compute_weights(min_root_x1[i], alpha)

# print('min_root_w1={}'.format(min_root_w1))
# print('min_root_x1={}'.format(min_root_x1))
# print('min_root_w2={}'.format(min_root_w2))
# print('min_root_x2={}'.format(min_root_alphas))

################################################################################

# Make plot of first root for 1/6 < alpha < alpha1. These are PO rules.
outside_alphas = np.concatenate((np.linspace(1./6 + 1.e-6, 0.17), np.linspace(0.17, alpha1_trig)))
outside_x1_first = np.zeros(len(outside_alphas))
outside_w1_first = np.zeros(len(outside_alphas))
outside_w2_first = np.zeros(len(outside_alphas))

for i in xrange(len(outside_alphas)):
    alpha = outside_alphas[i]
    roots = compute_x1_analytical(alpha)
    # print('alpha={}, roots={}'.format(alpha, roots))
    outside_x1_first[i] = roots[0].real
    outside_w1_first[i], outside_w2_first[i] = compute_weights(outside_x1_first[i], alpha)

# The weights are positive but the points are outside, so this is a "PO" rule
# print('outside_w1_first={}'.format(outside_w1_first))
# print('outside_x1_first={}'.format(outside_x1_first))
# print('outside_w2_first={}'.format(outside_w2_first))

# Make plot of second root for 1/6 < alpha < alpha1. These are PO rules.
outside_x1_second = np.zeros(len(outside_alphas))
outside_w1_second = np.zeros(len(outside_alphas))
outside_w2_second = np.zeros(len(outside_alphas))

for i in xrange(len(outside_alphas)):
    alpha = outside_alphas[i]
    roots = compute_x1_analytical(alpha)
    # print('alpha={}, roots={}'.format(alpha, roots))
    outside_x1_second[i] = roots[1].real
    outside_w1_second[i], outside_w2_second[i] = compute_weights(outside_x1_second[i], alpha)

# The weights are positive but the points are outside, so this is a "PO" rule
# print('outside_w1_second={}'.format(outside_w1_second))
# print('outside_x1_second={}'.format(outside_x1_second))
# print('outside_w2_second={}'.format(outside_w2_second))

################################################################################

# Note: for alpha1 < alpha < 1/5, the roots are imaginary, so we skip computing those.

# Compute roots for 1/5 < alpha < alpha2. These are NI rules.
alpha1_alpha2 = np.linspace(1./5, alpha2_trig)
# first root
alpha1_alpha2_x1_first = np.zeros(len(alpha1_alpha2))
alpha1_alpha2_w1_first = np.zeros(len(alpha1_alpha2))
alpha1_alpha2_w2_first = np.zeros(len(alpha1_alpha2))
# second root
alpha1_alpha2_x1_second = np.zeros(len(alpha1_alpha2))
alpha1_alpha2_w1_second = np.zeros(len(alpha1_alpha2))
alpha1_alpha2_w2_second = np.zeros(len(alpha1_alpha2))

for i in xrange(len(alpha1_alpha2)):
    alpha = alpha1_alpha2[i]
    roots = compute_x1_analytical(alpha)
    print('alpha={}, roots={}'.format(alpha, roots))
    # First root
    alpha1_alpha2_x1_first[i] = roots[0].real
    alpha1_alpha2_w1_first[i], alpha1_alpha2_w2_first[i] = compute_weights(alpha1_alpha2_x1_first[i], alpha)
    # Second root
    alpha1_alpha2_x1_second[i] = roots[1].real
    alpha1_alpha2_w1_second[i], alpha1_alpha2_w2_second[i] = compute_weights(alpha1_alpha2_x1_second[i], alpha)

# First root
# print('alpha1_alpha2_x1_first={}'.format(alpha1_alpha2_x1))
# print('alpha1_alpha2_w1_first={}'.format(alpha1_alpha2_w1_first))
# print('alpha1_alpha2_w2_first={}'.format(alpha1_alpha2_w2_first))

# Second root
# print('alpha1_alpha2_x1_second={}'.format(alpha1_alpha2_x1))
# print('alpha1_alpha2_w1_second={}'.format(alpha1_alpha2_w1_second))
# print('alpha1_alpha2_w2_second={}'.format(alpha1_alpha2_w2_second))

################################################################################

# Make plots

# Plot (alpha, x1(alpha)) for 0 < alpha < 1/6
fig = plt.figure()
ax1 = fig.add_subplot(111)
# Plot line x=1/6
ax1.plot([1./6,1./6], [0.44,0.55],color='lightgray', linestyle='--', linewidth=1)
# Plot line x=0
ax1.plot([0.,0.],[0.44,0.55], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=0.5
ax1.plot([-0.1,0.18], [0.5,0.5], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=min_x1
ax1.plot([-0.1,0.18], [min_x1,min_x1], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=r1
ax1.plot([-0.1,0.18], [r1,r1], color='lightgray', linestyle='--', linewidth=1)
ax1.plot(alphas, x1, color='black', marker=None)
ax1.plot([0], [r1], color='black', linestyle='', marker='o')
ax1.plot([1./6], [0.5], color='black', linestyle='', marker='o')
ax1.plot([min_alpha], [min_x1], color='black', linestyle='', marker='o')
ax1.plot([r2], [r1], color='black', linestyle='', marker='o')
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$x_1(\alpha)$')
ax1.text(0.+.002, r1+.001, 'PB:\,$(0,r_1)$')
ax1.text(1./6-.03, 0.5-.004, r'PB:\,$(\frac{1}{6}, \frac{1}{2})$')
ax1.text(r2, r1-.004, r'PI:\,$(r_2, r_1)$')
ax1.text(min_alpha-.013, min_x1+.0055, r'PI:')
ax1.text(min_alpha-.013, min_x1+.002, r'$(x_1(\alpha_3), \alpha_3)$')
ax1.set_xlim([-0.01, 0.175])
ax1.set_ylim([0.443, 0.505])
plt.savefig('plot_00020_x1_vs_alpha.pdf', format='pdf')

# Plot (x1(alpha), alpha) for 0 < alpha < 1/6
# We use the same values, simply reversing the x and y data.
# We make a separate plot so that the details of the first plot are not lost.
fig = plt.figure()
ax1 = fig.add_subplot(111)
# Plot line y=1/6
ax1.plot([0.4,0.51], [1./6,1./6], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=r2
# ax1.plot([0.4,0.51], [r2,r2], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=min_alpha
ax1.plot([0.4,0.51], [min_alpha,min_alpha], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=0
ax1.plot([0.4,0.51], [0,0], color='lightgray', linestyle='--', linewidth=1)
# Plot line x=min_x1
ax1.plot([min_x1,min_x1], [-1,1], color='lightgray', linestyle='--', linewidth=1)
# Plot line x=r1
ax1.plot([r1,r1], [-1,1], color='lightgray', linestyle='--', linewidth=1)
# Plot line x=0.5
ax1.plot([0.5,0.5], [-1,1], color='lightgray', linestyle='--', linewidth=1)
# Plot the data - note we just flip the x and y coordinates from the first plot,
# since the solutions are symmetric about the line y=x
ax1.plot(x1, alphas, color='black', marker=None)
ax1.plot([r1], [r2], color='black', linestyle='', marker='o')
ax1.plot([0.5], [1./6], color='black', linestyle='', marker='o')
ax1.plot([min_x1], [min_alpha], color='black', linestyle='', marker='o')
ax1.plot([r1], [0], color='black', linestyle='', marker='o')
ax1.text(r1+.001, r2-.008, r'PI:\,$(r_1, r_2)$')
ax1.text(min_x1+.001, min_alpha+.0015, r'PI:\,$(\alpha_3, x_1(\alpha_3))$')
ax1.text(0.5 - .005, 1./6 - .015, r'PB:\,$(\frac{1}{2}, \frac{1}{6})$')
ax1.text(r1+.001, 0+.001, r'PB:\,$(r_1, 0)$')
ax1.set_xlim([min_x1 - .005, 0.505])
ax1.set_ylim([-0.02, 0.18])
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$x_1(\alpha)$')
plt.savefig('plot_00020_x1_vs_alpha_2.pdf', format='pdf')

### # Plot w1 and w2 vs. alpha. Do they cross somewhere?
### fig = plt.figure()
### ax1 = fig.add_subplot(111)
### # Plot line y=1/12. The weights are symmetric about this line.
### ax1.plot([-0.1,0.18], [1./12,1./12], color='lightgray', linestyle='--', linewidth=1)
### # Plot vertical line for the value of alpha where the weights are equal.
### ax1.plot([alpha_weights_equal,alpha_weights_equal], [0.,0.2],color='lightgray', linestyle='--', linewidth=1)
### ax1.plot(alphas, w1, color='blue', marker=None, linewidth=2, label=r'$w_1$')
### ax1.plot(alphas, w2, color='red', marker=None, linewidth=2, label=r'$w_2$')
### ax1.set_xlabel(r'$\alpha$')
### ax1.set_xlim([-0.01, 0.175])
### ax1.set_ylim([0.01, 0.155])
### ax1.legend()
### plt.savefig('plot_00020_weights_vs_alpha.pdf', format='pdf')

# Plot "min" root NI/NB branch (alpha, x1(alpha)) for 0 < alpha < min_alpha
# and then the symmetric part, since it is "connected" in this case.
fig = plt.figure()
ax1 = fig.add_subplot(111)
# Plot dashed line for y=x
# ax1.plot([0,1], [0,1], color='lightgray', linestyle='--', linewidth=1)
# Plot dashed line for y=0
# ax1.plot([0,1], [0,0], color='lightgray', linestyle='--', linewidth=1)
# Plot dashed line for x=0
# ax1.plot([0,0], [0,1], color='lightgray', linestyle='--', linewidth=1)
# line y=x1(alpha3)
ax1.plot([0,1], [x1_alpha3,x1_alpha3],color='lightgray', linestyle='--', linewidth=1)
# line x=x1(alpha3)
ax1.plot([x1_alpha3,x1_alpha3], [0,1], color='lightgray', linestyle='--', linewidth=1)
# Plot the data, which forms a symmetric path about the origin.
# We plot the first half of data up to the singular point, then
# swap the arrays to plot the other half.
ax1.plot(min_root_alphas, min_root_x1, color='black', linestyle='--', marker=None)
ax1.plot(min_root_x1, min_root_alphas, color='black', linestyle='--', marker=None)
# Plot single points
ax1.plot([0], [r2], color='black', linestyle='', marker='o')
ax1.plot([r2], [0], color='black', linestyle='', marker='o')
ax1.plot([min_alpha], [min_alpha], color='black', linestyle='', marker='o', markerfacecolor='white')
# Add plot labels
ax1.text(0+.002, r2+.001, r'NB:\,$(0, r_2)$')
ax1.text(r2-.035, 0+.001, r'NB:\,$(r_2, 0)$')
# ax1.text(min_alpha-.05, min_alpha-.005, r'$(x_1(\alpha_3), x_1(\alpha_3))$')
ax1.text(0.01, x1_alpha3+.002, r'$x_1(\alpha_3)$')
ax1.axis('square')
ax1.set_xlim([-0.001, r2 + .01])
ax1.set_ylim([-0.001, r2 + .01])
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$x_1(\alpha)$')
plt.savefig('plot_00020_min_root.pdf', format='pdf')

### # Plot weights on "min" root NI/NB branch (w1 vs. alpha, w2 vs. alpha)
### fig = plt.figure()
### ax1 = fig.add_subplot(111)
### # Plot line y=1/12. The weights are symmetric about this line.
### # ax1.plot([-0.1,0.18], [1./12,1./12], color='lightgray', linestyle='--', linewidth=1)
### # Plot vertical line for the value of alpha where the weights are equal.
### # ax1.plot([alpha_weights_equal,alpha_weights_equal], [0.,0.2],color='lightgray', linestyle='--', linewidth=1)
### ax1.plot(min_root_alphas, min_root_w1, color='blue', marker=None, linewidth=2, label=r'$w_1$')
### ax1.plot(min_root_alphas, min_root_w2, color='red', marker=None, linewidth=2, label=r'$w_2$')
### # Plot second half of data, which is just the other weight due to symmetry, and
### # we can use the min_root_x1 for the x-values for the same reason.
### ax1.plot(min_root_x1, min_root_w2, color='blue', marker=None, linewidth=2)
### ax1.plot(min_root_x1, min_root_w1, color='red', marker=None, linewidth=2)
### ax1.set_xlabel(r'$\alpha$')
### # ax1.set_xlim([-0.01, 0.175])
### # ax1.set_ylim([0.01, 0.155])
### ax1.legend()
### plt.savefig('plot_00020_min_root_weights_vs_alpha.pdf', format='pdf')

# Make plot for 1/6 < alpha < alpha1_trig. For these cases there are
# two solution branches, both of which are outside the reference element.
# This is a PO rule, so we use a solid line.
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(outside_alphas, outside_x1_first, color='black', linestyle='-', marker=None)
ax1.plot(outside_alphas, outside_x1_second, color='black', linestyle='-', marker=None)
# line x=alpha_1
ax1.plot([alpha1_trig,alpha1_trig], [0,3],color='lightgray', linestyle='--', linewidth=1)
# line x=1/6
ax1.plot([1./6,1./6], [0,3],color='lightgray', linestyle='--', linewidth=1)
# line y=0.5
ax1.plot([0,1], [0.5,0.5],color='lightgray', linestyle='--', linewidth=1)
# Plot single points
ax1.plot([1./6], [0.5], color='black', linestyle='', marker='o')
ax1.plot([alpha1_trig], [x1_alpha1], color='black', linestyle='', marker='o')
# Point labels
ax1.text(1./6+5.e-5, 0.5+.05, r'PB:\,$(\frac{1}{6},\frac{1}{2})$')
ax1.text(alpha1_trig-.0011, x1_alpha1, r'PO:\,$(\alpha_1, x_1(\alpha_1))$')
# Labels, limits, and legends
ax1.set_xlim([1./6-.0001, 0.171])
ax1.set_ylim([0.4, 2])
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$x_1(\alpha)$')
plt.savefig('plot_00020_outside_alphas.pdf', format='pdf')
