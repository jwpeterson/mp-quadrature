function [q,w] = gquad(n)

% [q,w] = gquad(n)
% Computes the points and weights of
% an n-point Guassian quadrature rule.
% It uses the function leg(n) to compute
% the Legendre polynomials up to order
% to determine the quadrature points.

% First get the roots of the desired Legendre function.
[P,R] = leg(n);
q = R{length(R)};



% Now compute the n weights
for j=1:n
  % Form the interpolating Lagrange polynomial
  % by successively convolving and the normalization
  % constant by multiplication.
  lagrange = 1;
  lagconst = 1;
  for i=1:n
    if i ~= j
      lagrange = conv([1 -q(i)], lagrange);
      lagconst = lagconst * (q(j) - q(i));
    end
  end
 %lagrange
 %lagconst

 % Compute the integral of the Lagrange polynomial
 % between -1 and 1, and scale by the constant to
 % get the weight.
 intlagrange = polyint(lagrange);
 w(j) = 1/lagconst * (polyval(intlagrange, 1) - ...
		      polyval(intlagrange,-1));
end

msg = ['The sum of the weights is: ', num2str(sum(w)), '.'];
disp(msg);
w=w';
