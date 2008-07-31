function [x,w] = jacobi_rule(alpha,beta,n)
% This routine computes the points and weights for
% a Jacobi quadrature rule with n points.  It relies
% on the jacobi_points routine to first determine the
% roots of P_n^(a,b)(x).  Then it computes the weights
% according to the formula given by Stroud.  This allows
% the point and weight finding routines to be done
% separately, but we do need to compute the polynomial
% values one extra time.

% Set number printing format
format long e

% First get the points
x = jacobi_points(alpha, beta, n);

% Get the derivative dP_n and the (n-1)^st function value at
% all of the quadrature points.
[p,dp,pnm1] = jacobi_recur(alpha, beta, n, x);

% Get the jacobi coefficients for j=1:n.
% We need the c's to compute the weight formula.
[b,c] = jacobi_constants(alpha, beta, linspace(1,n,n));

% For the "true" sum of the weights, Stroud gives (I need to look this up...):
tsa = 2^(alpha+beta+1) * factorial(alpha) / factorial(alpha+beta+1) * factorial(beta);

% Compute the common weight coefficient.  This multiplies each weight.
% We do it with a loop to avoid problems with overflow.  We can
% do this b/c we assume alpha and beta are integers.  This coefficient
% is = tsa*prod(c(2:n)).  Since we can scale the weights to sum up to
% anything we want, is multiplying by this extra tsa value really
% necessary?
%% loop_max = max(alpha+beta+1, n);
%% weight=1.;
%% for j=1:loop_max
%%   if (j<=alpha)
%%     weight = weight * j; % alpha!
%%   end
%% 
%%   if (j<=alpha+beta+1)
%%     weight = weight * (2. / j); % 2^(alpha+beta+1) / (alpha+beta+1)!
%%   end
%%   
%%   if (j<=beta)
%%     weight = weight * j; % beta!
%%   end
%% 
%%   if ((j>=2) & (j<=n))
%%     weight = weight * c(j); % c(2)*c(3)*...*c(n)
%%   end
%% end

%weight

% Finally, to get the individual weights, we divide the weight coefficient by
% P_{n-1}(x_i) P'_n(x_i)
%w = weight ./ dp ./ pnm1;
w = prod( c(2:n) ) ./ dp ./ pnm1;

% Also print out the sum of weights before returning
disp(['Sum of weights=', num2str(sum(w), '%16.15e')]);
