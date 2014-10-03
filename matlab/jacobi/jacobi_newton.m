function r = jacobi_newton(a, b, n, x)
% Uses Newton's method to converge to a root of the
% Jacobi polynomial P(a,b,n) given initial guess x.
%
% * Newton is convenient because the derivative
%   is readily available from the recurrence relation.
% * This routine, like jacobi_recur(), expects its
%   x-arguments to lie in [-1 1]

% Tolerance to converge Newton residual to.
tol = 1.e-15;

% Current residual of the Newton system.
residual = 1.;
n_its=0;
max_its=30;

while ((residual > tol) & (n_its<max_its))
  [p, dp, pnm1] = jacobi_recur(a,b,n,x);
  x = x - p ./ dp;
  residual = norm(p);
  n_its=n_its+1;
end

if ((n_its==max_its) & (residual > tol))
  error('Iterations did not converge!');
end

r=x;
