import numpy as np
import matplotlib.pyplot as plt
import sys
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
Given alpha, computes x1(alpha), throwing an error if an invalid
value is obtained.
"""
def compute_x1(alpha):
    # Don't accept input values outside the range
    if (alpha <= 0.) or (alpha >= 0.5):
        print('Invalid input, 0 < alpha < 0.5 is required.')
        sys.exit(1)

    sigma = f(alpha) / g(alpha)
    roots = np.roots([15*sigma, -(9*sigma + 6), (sigma + 4), -0.5])

    # We want to pick the positive, real root which is in the interval
    # [0, 1/2] and which is farthest from the input value of alpha. I didn't
    # realize this before, but alpha itself should always be one of the roots?
    # print('alpha={}'.format(alpha))
    # print('roots={}'.format(roots))

    # Return value
    x1 = 0
    found_x1 = False
    for candidate_root in roots:
        if np.isreal(candidate_root) and (np.abs(candidate_root - alpha) > 1.e-3) \
           and (candidate_root > 0.) and (candidate_root < 0.5):
            found_x1 = True
            x1 = candidate_root
            break

    # Error if an acceptable root was not found.
    if not found_x1:
        print('Invalid root = {} found for alpha={}.'.format(x1, alpha))
        sys.exit(1)

    return x1

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
    x1 = compute_x1(alpha)
    w1, w2 = compute_weights(x1, alpha)
    return w1 - 1./12

################################################################################

# Some useful constants
r1 = (9 + np.sqrt(21))/30
r2 = (9 - np.sqrt(21))/30

################################################################################

# There should be a "symmetric" part for larger values of alpha > (9 +
# sqrt(21))/30 as well. When searching for the correct root, we want
# the one that is furthest from the input value alpha, not closest
# to any fixed value.
# alpha_vec = np.linspace(r1 + 1.e-6, .5 - 1.e-6)
# for alpha in alpha_vec:
#     x1 = compute_x1(alpha)
#     w1, w2 = compute_weights(x1, alpha)
#
#     # Print current result
#     print('---')
#     print('w1={}'.format(w1))
#     print('x1={}'.format(x1))
#     print('w2={}'.format(w2))
#     print('x2={}'.format(alpha))
#
# # Early
# sys.exit(0)

################################################################################

# Find the minimum value of x1(alpha). Interestingly, it does not
# seem to coincide with either of the limits alpha=0 or alpha=1/6.
# method='L-BFGS-B'
# method='Newton-CG' # requires Jacobian
# method='CG' # cannot handle constraints or bounds
result = minimize(compute_x1, 0.11, method='CG', \
                  options={'disp': False,
                           'gtol' : 1.e-10})

# Extract the min x1 value and the alpha where it occurs.
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
    x1_weights_equal = compute_x1(alpha_weights_equal)
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

# The curve goes from:
# (x1,x2) = ((9 + np.sqrt(21))/30, 0) to
# (x1,x2) = (1/2, 1/6)
# But w1/w2 -> \infty at alpha=1/6
alphas1 = np.linspace(1.e-6, 0.15)
alphas2 = np.linspace(0.15, 1./6 - 1.e-6)
alphas = np.concatenate((alphas1,alphas2))
# print('alphas={}'.format(alphas))
w1 = np.zeros(len(alphas))
w2 = np.zeros(len(alphas))
x1 = np.zeros(len(alphas))

for i in xrange(len(alphas)):
    alpha = alphas[i]
    x1[i] = compute_x1(alpha)
    w1[i], w2[i] = compute_weights(x1[i], alpha)

# Print results
# print('w1={}, w2={}'.format(w1, w2))

# Plot results
fig = plt.figure()
ax1 = fig.add_subplot(111)

# Plot y=x line of symmetry
ax1.plot([0,.2],[0,.2], color='lightgray', linestyle='--', linewidth=1)

# Plot line x=0
ax1.plot([0,0],[0,.2], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=0
ax1.plot([0,.2],[0,0], color='lightgray', linestyle='--', linewidth=1)

# Plot main solutions (and symmetric solutions
ax1.plot(w1, w2, color='black', marker=None)
# ax1.plot(w2, w1, color='salmon', marker=None)

# Plot other solutions.

# The PI solution: should be on the curve since we are expanding about
# the PI solution.
ax1.plot([1./12 - np.sqrt(21)/168],
         [1./12 + np.sqrt(21)/168], color='red', linestyle='', marker='o')

# The PB/NB solution pair
ax1.plot([(39 - np.sqrt(21))/240],
         [(1 + np.sqrt(21))/240], color='blue', linestyle='', marker='o')
ax1.plot([(1 - np.sqrt(21))/240],
         [(39 + np.sqrt(21))/240], color='green', linestyle='', marker='o')

# The PB solution without partner
ax1.plot([1./60],
         [3./20], color='khaki', linestyle='', marker='o')

#ax1.set_xlim([0, 0.1])
#ax1.set_ylim([0, 0.1])
ax1.set_xlabel(r'$w_1$')
ax1.set_ylabel(r'$w_2$')
ax1.axis('equal')
plt.savefig('plot_00020.pdf', format='pdf')

# Plot (alpha, x1(alpha))
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
ax1.plot(alphas, x1, color='black', marker=None)
ax1.plot([0], [r1], color='black', linestyle='', marker='o')
ax1.plot([1./6], [0.5], color='black', linestyle='', marker='o')
ax1.plot([min_alpha], [min_x1], color='black', linestyle='', marker='o')
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$x_1(\alpha)$')
ax1.text(0.+.002, r1+.001, r'$\alpha = 0$')
ax1.text(0.+.002, r1+.0035, 'PB,')
ax1.text(1./6-.018, 0.5-.004, r'$\alpha = \frac{1}{6}$')
# ax1.text(1./6-.018, 0.5-.001, 'PB')
ax1.text(1./6-.03, 0.5-.00325, 'PB,')
ax1.text(min_alpha-.01, min_x1+.002, r'$\alpha \approx 0.109$')
ax1.set_xlim([-0.01, 0.175])
ax1.set_ylim([0.443, 0.505])
plt.savefig('plot_00020_x1_vs_alpha.pdf', format='pdf')

# Plot w1 and w2 vs. alpha. Do they cross somewhere?
fig = plt.figure()
ax1 = fig.add_subplot(111)
# Plot line y=1/12. The weights are symmetric about this line.
ax1.plot([-0.1,0.18], [1./12,1./12], color='lightgray', linestyle='--', linewidth=1)
# Plot vertical line for the value of alpha where the weights are equal.
ax1.plot([alpha_weights_equal,alpha_weights_equal], [0.,0.2],color='lightgray', linestyle='--', linewidth=1)
ax1.plot(alphas, w1, color='blue', marker=None, linewidth=2, label=r'$w_1$')
ax1.plot(alphas, w2, color='red', marker=None, linewidth=2, label=r'$w_2$')
ax1.set_xlabel(r'$\alpha$')
ax1.set_xlim([-0.01, 0.175])
ax1.set_ylim([0.01, 0.155])
ax1.legend()
plt.savefig('plot_00020_weights_vs_alpha.pdf', format='pdf')
