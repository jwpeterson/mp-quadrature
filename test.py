#!/usr/bin/env python

import numpy as np
from scipy.optimize import fsolve

# This is a hard-coded solver for the degree=3 Ro3-invariant case with
# two median orbits. This case behaved very strangely in our solver,
# seeming to have many, many solutions.
def residual(x):
    r = np.zeros(4)
    for i in xrange(2):
        wi = x[2*i]
        xi = x[2*i+1]
        # (0,0)
        r[0] += 3*wi
        # (2,0)
        r[1] += wi * (6*xi**2 - 4*xi + 1)
        # (3,0)
        r[2] += wi * (-6*xi**3 + 12*xi**2 - 6*xi + 1)
        # (2,1)
        r[3] += wi * xi * (3*xi**2 - 3*xi + 1)

    r = np.subtract(r, [1./2, 1./12, 1./20, 1./60])
    return r

# Computes the Jacobian analytically
def jacobian(x):
    J = np.zeros((4,4))

    # First row is constant
    J[0,0] = 3
    J[0,2] = 3

    for i in xrange(2):
        wi = x[2*i]
        xi = x[2*i+1]

        # Row 1
        J[1, 2*i] = 6*xi**2 - 4*xi + 1 # deriv wrt wi
        J[1, 2*i+1] = wi * (12*xi - 4) # deriv wrt xi

        # Row 2
        J[2, 2*i] = -6*xi**3 + 12*xi**2 - 6*xi + 1 # deriv wrt wi
        J[2, 2*i+1] = wi * (-18*xi**2 + 24*xi - 6) # deriv wrt xi

        # Row 3
        J[3, 2*i] = xi * (3*xi**2 - 3*xi + 1) # deriv wrt wi
        J[3, 2*i+1] = wi * (27*xi**2 - 6*xi + 1) # deriv wrt xi

    return J

# Solution vector is in the order: [w1, x1, w2, x2]

# Random initial guess
lower_limits = np.zeros(4);
upper_limits = np.array([1./6, 0.5, 1./6, 0.5])
x = np.random.uniform(low=lower_limits, high=upper_limits)

# Some random initial guesses that converged
# x = np.array([0.08674254, 0.20938463, 0.14345518, 0.42824847])
# x = np.array([0.08691546, 0.21344877, 0.15825936, 0.08199257]) # converges to a negative wt soln
# x = np.array([0.14784719, 0.29403371, 0.06712839, 0.00667038]) # converges to a point outside

# What the first random initial guesses above converged to
# x = np.array([0.1124398, 0.14837781, 0.05422687, 0.4535146])

# Some solutions found by the C++ code (for verification)
# Note: the first initial guess doesn't converge to itself, but to
# some very nearby solution, which is also incredibly weird to me.
# x = np.array([3.538027752991e-2, 4.516639556788e-2, 1.314470832136e-1, 4.462746877117e-1])
# x = np.array([5.543043859271e-2, 9.232121925129e-2, 1.112362280739e-1, 4.459138996659e-1])
# x = np.array([5.711759998754e-2, 4.523356909977e-1, 1.095490666791e-1, 1.465791804406e-1])
# x = np.array([7.793391579768e-2, 1.212261086016e-1, 8.873275086898e-2, 4.458533573224e-1])
# x = np.array([1.204263852130e-1, 4.468905352570e-1, 4.624028145360e-2, 7.538511462816e-2])
# x = np.array([1.245621175415e-1, 1.552459464812e-1, 4.210454912515e-2, 4.603355008546e-1])
# x = np.array([7.526226664663e-2, 1.184342266352e-1, 9.140440002003e-2, 4.456902313190e-1])
# x = np.array([1.333490430034e-1, 4.494216822955e-1, 3.331762366318e-2, 4.175729528808e-2])
# x = np.array([1.485869152981e-1, 1.661360798447e-1, 1.807975136854e-2, 4.955736994366e-1])
# x = np.array([1.384772362022e-1, 4.509265487983e-1, 2.818943046438e-2, 2.300511903023e-2])

# And there are many, many more!
# x = np.array([
# 8.7968789711389793435006162846315e-2,
# 4.4590889747817179822420514240849e-1,
# 7.8697876955276873231660503820351e-2,
# 1.2200025851717823099515109050002e-1])

# x = np.array([
# 4.0530576752184521725793976191672e-2,
# 6.2344546547003829041672350905925e-2,
# 1.2613608991448214494087269047499e-1,
# 4.4781928344493530500276639517784e-1])

# x = np.array([
# 2.5771715771331835575174791068785e-2,
# 1.2468756236904283744035960268095e-2,
# 1.4089495089533483109149187559788e-1,
# 4.5176875855127644088192186961460e-1])

# x = np.array([
# 2.1585379122821437987049726854104e-2,
# 4.8658154284247656262552542413326e-1,
# 1.4508128754384522867961693981256e-1,
# 1.6476150682211463282117998527133e-1])

# x = np.array([
# 8.6677262602444517186824380088901e-2,
# 4.4601196987393666566496242143040e-1,
# 7.9989404064222149479842286577766e-2,
# 1.2328567763035714822089586239362e-1])

# x = np.array([
# 1.3514500636744604509230060616125e-1,
# 4.4990940050485156545552288382557e-1,
# 3.1521660299220621574366060505420e-2,
# 3.5665544866737465105449499598089e-2])

# x = np.array([
# 4.9793911745574969085127148954628e-2,
# 4.5561683367407423165179141931074e-1,
# 1.1687275492109169758153951771204e-1,
# 1.5101027921916621385589628328530e-1])

# x = np.array([
# 3.1386024253318503689859605439649e-2,
# 3.5186328882414804056159504067558e-2,
# 1.3528064241334816297680706122702e-1,
# 4.4994786568640857034903149408922e-1])

# x = np.array([
# 1.8678342730319004358105067898138e-2,
# 4.9385495412351989975299305916744e-1,
# 1.4798832393634766230856159876853e-1,
# 1.6590709891340603950320438271269e-1])

# x = np.array([
# 5.1905167571475480183119331518242e-2,
# 8.6312483542116624941053333998057e-2,
# 1.1476149909519118648354733514842e-1,
# 4.4621849457434170156811154555106e-1])

# If I use the finite differenced Jacobian option (default) in fsolve,
# it tells me the Jacobian is nearly singular at root, i.e. det(J) = -3.90337195411e-17
# If I use my analytically computed Jacobian, it is not singular. Both of them
# converge to the same answer. What is going on?
x = np.array([8.796878971138e-2, 4.459088974781e-1, 7.869787695527e-2, 1.220002585171e-1])

print('initial guess={}'.format(x))

r = residual(x)
print('initial residual={}'.format(r))

# If iflag==1 if solution found, otherwise, check contents of "mesg"
(sol, infodict, iflag, mesg) = fsolve(residual, x, fprime=jacobian, full_output=True)

if iflag == 1:
    # The solution (if converged) is in the 0'th entry of the returned tuple
    print('Solution = {}'.format(sol))
    print('Number of function evaluations = {}'.format(infodict['nfev']))
    print('Final residual = {}'.format(infodict['fvec']))
    print('Final residual norm = {}'.format(np.linalg.norm(infodict['fvec'])))

    Q = np.array(infodict['fjac'])
    Rvec = infodict['r']
    # print('Upper-triangular R from QR of final Jacobian = \n{}'.format(Rvec))

    # Convert from upper triangular Rvec -> into 4x4 matrix, R.
    # The function triu_indices returns the triangular indices of an upper-triangular matrix.
    R = np.zeros((4,4))
    (rows, cols) = np.triu_indices(4)
    for k in xrange(len(rows)):
        R[rows[k], cols[k]] = Rvec[k]

    # print('Orthogonal Q from QR of final Jacobian = \n{}'.format(Q))
    # print('Upper-triangular R from QR of final Jacobian = \n{}'.format(R))

    # Compute the Schur form of J by multiplying Q and R together. The
    # diagonal entries are the eigenvalues of the original matrix.
    JacQR = np.multiply(Q, R)
    # print('QR = Final Jacobian in Schur form = \n{}'.format(JacQR))

    # Extract the diagonal of R
    diagR = np.zeros(4)
    for k in xrange(4):
        diagR[k] = R[k,k]

    # The determinant of the Jacobian is therefore product of the
    # diagonal entries of R.
    det = np.product(diagR)
    # Condition number of J is abs(lambda_max) / abs(lambda_min)
    condJ = np.max(np.abs(diagR)) / np.min(np.abs(diagR))

    print('det(QR of J) = {}'.format(det))
    print('cond(QR of J) = {:.10e}'.format(condJ))

    # Compute J analytically at the solution
    Jac = jacobian(sol)
    print('Final Jacobian (analytical) = \n{}'.format(Jac))

    # Compute determinant of analytical Jacobian
    det = np.linalg.det(Jac)
    condJ = np.linalg.cond(Jac)
    (lam, v) = np.linalg.eig(Jac)
    print('det(J) = {}'.format(det))
    print('cond(J) = {:.10e}'.format(condJ))
    print('eig(J) = {}'.format(lam))

    # Wolfram-alpha input for transformed variables (eta1, eta2, eta3)
    # row reduce {{-4, 6, 0, -1/12}, {-6, 12, -6, -7/60}, {1, -3, 3, 1/60}}
    # Wolfram-alpha input for the next larger problem (degree 5, 7 unknowns)
    # This seemed to be senstive to the amount of space and the newlines (?) so
    # I just removed all of them.
    # row reduce {{-2/9,-4,6,0,0,0,-1/12},{-8/27,-6,12,-6,0,0,-7/60},{1/27,1,-3,3,0,0,1/60},{-26/81,-8,24,-32,18,0,-2/15},{-80/243,-10,40,-80,80,-30,-1/7},{1/243,1,-8,24,-31,15,1/210}}
    # For the (1,1,0,1,0) configuration with 7 QPs. I'm not sure if this is approach really makes sense, but Wolfram-alpha tells us this matrix is full rank (4).
    # row reduce {{1,3,3,0,0,0,1/2},{1/9,1,1,-4,6,0,1/12},{1/27,1,1,-6,12,-6,1/20},{1/27,1,0,1,-3,3,1/60}}

    # Compute the transformed variables according to their definitions.
    eta3 = sol[0]*sol[1]**3 + sol[2]*sol[3]**3
    eta2 = sol[0]*sol[1]**2 + sol[2]*sol[3]**2
    eta1 = sol[0]*sol[1] + sol[2]*sol[3]
    print('eta1={}'.format(eta1))
    print('eta2={}'.format(eta2))
    print('eta3={}'.format(eta3))

    # These differences should be approximately zero if our linear algebra was correct!
    print('eta1 - (1./40 + 3 * eta3)={}'.format(eta1 - (1./40 + 3 * eta3)))
    print('eta2 - (1./360 + 2 * eta3)={}'.format(eta2 - (1./360 + 2 * eta3)))
