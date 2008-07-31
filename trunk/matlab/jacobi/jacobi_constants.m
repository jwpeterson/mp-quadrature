function [b,c] = jacobi_constants(alpha,beta,n)
% This routine computes and returns the recurrence
% constants b_n and c_n for the Jacobi polynomials given
% by Stroud.  n may be a vector of values, in which case
% b and c will be returned as vectors as well.

% I tried interleaving multiplications and divisions but I'm not sure it matters
b = (alpha+beta) ./ (alpha+beta+2.*n) .* (beta-alpha) ./ (alpha+beta+2.*n-2.);
c = 4.*(n-1.) ./ (alpha+beta+2.*n-1.) .* (alpha+n-1.) ./ (alpha+beta+2.*n-2.).^2 .* (beta+n-1.) ./ (alpha+beta+2.*n-3.) .* (alpha+beta+n-1.);
