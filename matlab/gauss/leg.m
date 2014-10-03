function [P,R] = leg(n)

% [P,R] = leg(n)
% Computes Pn, the Legendre polynomial
% of order n on the interval [-1,1].
%
% Returns:
% P -- the coefficients of
%      the Legendre polynomials up to order n
% R -- the roots of the corresponding Legendre
%      polynomial



% The recursion relation needs the first
% two polynomials to get started.
P{1} = 1;
P{2} = [1 0];

R{1} = [];
R{2} = 0;

% If the user wants the constant or linear
% Legendre polynomials, we can return early.
if (n <= 1)
  return
end


% Loop to compute the polynomials.  When
% accessing the array P, make sure to increment
% the index by 1 every time.
for l=1:n-1
  term1 = conv([2*l+1 0], P{l+1});
  term2 = -l*P{l-1+1};
  term3 = l+1;

  % The lengths of term1 and term2 must match
  % term1 will always be longer.  Therefore we
  % pad term2 with zeros.
  term2 = [zeros(1,length(term1) - length(term2)), term2];

  % Compute the next polynomial
  P{l+1+1} = (term1 + term2) / term3;
end

for i=3:length(P)
  R{i} = sort(roots(P{i}));
end
