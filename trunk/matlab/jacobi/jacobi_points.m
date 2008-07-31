function r = jacobi_points(a, b, n)
% This routine implements Stroud's initial guesses for computing the
% roots of the Jacobi polynomial P(a,b,n).  The initial guesses are
% combined with Newton's method (see jacobi_newton.m) which iterates
% until the desired tolerance is reached.
%
% Note: the initial guess for each root is
% *based on the Newton solution to the previous root*
% This means our Newton call can't vectorize, but that's OK...
  
% The result vector.  The roots will be returned in this.  
r = zeros(n,1); 

% Some variable combinations we will need frequently...
an = a / n;
bn = b / n;

% The roots are computed in reverse order starting with the root
% closest to x=1.  The initial guess for each subsequent root
% depends on the previous converged solution.
for j=1:n
  
  % Stroud's "largest zero," the root closest to x=1
  if (j==1)
    R1 = (1.+a)*(2.78/(4.+n*n) + .768*an/n);
    R2 = 1. + 1.48*an + .96*bn + .452*an*an + .83*an*bn;
    r(1) = 1.- R1/R2;
  end
  
  % Stroud's "second zero"
  % Use converged r(1) to generate a guess for the second root, r(2)
  if (j==2)
    R1 = (4.1+a)/((1.+a)*(1.+.156*a));
    R2 = 1. + .06*(n-8.)*(1.+.12*a)/n;
    R3 = 1. + .012*b*(1. + .25*abs(a))/n;
    r(2) = r(1) - R1*R2*R3*(1.-r(1));
  end
  
  % Stroud's "third zero"
  % Use converged r(2) to generate a guess for r(3) and r(2)
  if (j==3)
    R1 = (1.67 + .28*a)/(1. + .37*a);
    R2 = 1. + .22*(n-8.)/n;
    R3 = 1. + 8.*b/((6.28 + b)*n*n);
    r(3) = r(2) - R1*R2*R3*(r(1) - r(2));
  end
  
  % Stroud's "middle zeros"
  if ((j > 3) & (j<n-1))
    r(j) = 3.*r(j-1) - 3.*r(j-2) + r(j-3);
  end
  
  % Stroud's "second last zero"
  if ((j==n-1) & (n>3))
    R1 = (1. + .235*b)/(.766 + .119*b);
    R2 = 1. / (1. + .639*(n-4.)/(1. + .71*(n-4.)));  
    R3 = 1. / (1. + 20.*a/((7.5+a)*n*n) );
    r(j) = r(j-1) + R1*R2*R3*(r(j-1) - r(j-2));
  end
  
  % Stroud's "last zero"
  if ((j==n) & (n>3))
    R1 = (1. + .37*b)/(1.67 + .28*b);
    R2 = 1. / (1. + .22*(n-8.)/n);
    R3 = 1. / (1. + 8.*a/((6.28+a)*n*n));
    r(j) = r(j-1) + R1*R2*R3*(r(j-1) - r(j-2));
  end

  % Call the newton routine using the most-recently computed
  % initial guess.
  r(j) = jacobi_newton(a, b, n, r(j));

end %for j=1:n
  
% Put roots in order from smallest to largest
r=sort(r);
