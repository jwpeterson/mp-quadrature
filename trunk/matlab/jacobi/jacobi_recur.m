function [p,dp,pnm1] = jacobi_recur(a, b, n, x)
% This function computes the values of the Jacobi
% polynomial P(alpha,beta,n) and its derivative at x
% using the recursion given by Stroud:
%
% P_n = (x-b_n)*P_{n-1} - c_n*P_{n-2} , n >= 2
%
% where b_n = (a+b)(b-a) / (a+b+2n) / (a+b+2n-2)
% where c_n = 4*(n-1)*(a+n-1)*(b+n-1)*(a+b+n-1) / (a+b+2*n-1) / (a+b+2*n-2)^2 / (a+b+2*n-3)
%
% The polynomials are defined on [-1 1].  The routine is initialized by
% P_0 = 1
% P_1 = x + (a-b)/(a+b+2)
%
% And the coefficients b_n and c_n take care of the rest of the scaling.
% * Uses for this function include: Newton's Method for root finding.
% * It may be more accurate than pre-computing all n+1 polynomial coeffs
%   and evaluating the polynomial that way.
% * This works with x scalar or an array.
% See: Stroud and Secrest, "Gaussian Quadrature Formulas," 1966, Prentice Hall
%
% The Jacobi polynomials are typically normalized such that:
% P_n^(a,b)(1) = (n+a)! / n! / a!
% I don't think Stroud's recurrence relation follows this convention, but
% that scaling does not change the roots of course.
  
% Special return for n=0 and 1
if (n==0)
  p=1.;
  dp=0.;
  pnm1=0.; % ? This is just to avoid errors about unassigned vars.
  return;
end

if (n==1)
  p=0.5*(a+b+2)*x + 0.5*(a-b);
  dp=0.5*(a+b+2);
  pnm1=1.;
  return;
end


% Initialize for Stroud's recurrence relation

% P_{n-2} and P'_{n-2}
pnm2  = 1.;
dpnm2 = 0.;

% P_{n-1} and P'_{n-1}
pnm1  = x + (a-b)/(a+b+2.);
dpnm1 = 1.;



for j=2:n
  % Compute the recurrence coefficients b_j, c_j.
  [bj,cj] = jacobi_constants(a,b,j);
  
  % Compute P_{n}, dP/dx_{n} using the recurrence relation
  p  = (x-bj).*pnm1 - cj*pnm2;
  dp = (x-bj).*dpnm1 + pnm1 - cj*dpnm2;

  % At the last step, return before updating the values
  if (j==n)
    return;
  end
  
  % Update polynomial values
  pnm2 = pnm1;
  pnm1 = p;
  
  % Update polynomial derivs
  dpnm2 = dpnm1;
  dpnm1 = dp;
end
