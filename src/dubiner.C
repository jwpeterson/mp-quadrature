#include <cstdlib> // std::abort
#include "dubiner.h"
#include "common_definitions.h"
#include "conical.h"

void Dubiner::p_numeric(unsigned d,
                        const mpfr_class & xi,
                        const mpfr_class & eta,
                        std::vector<mpfr_class> & vals)
{
  // Make sure there are no old values in the vector.
  vals.clear();

  for (unsigned i=0; i<=d; ++i)
    for (unsigned j=0; j<=d; ++j)
      {
        // If degree is too high, continue
        if (i+j > d)
          continue;

        // Map (xi, eta) to barycentric coordinates
        mpfr_class
          zeta0 = xi,
          zeta1 = eta,
          zeta2 = 1 - xi - eta;

        // Compute transformed coordinate for evaluating P_i
        mpfr_class transformed1 = (zeta1-zeta0) / (zeta1+zeta0);

        // Compute P_i^(0,0)
        mpfr_class P_i = this->jacobi(/*n=*/i, /*alpha=*/0, /*beta=*/0, /*x=*/transformed1);

        // Scale P_i appropriately
        P_i *= pow(zeta1 + zeta0, i);

        // Compute transformed coordinate for evaluating P_j
        mpfr_class transformed2 = 2.*zeta2 - 1.;

        // Compute P_j^(2*i+1,0)
        mpfr_class P_j = this->jacobi(/*n=*/j, /*alpha=*/2*i+1, /*beta=*/0, /*x=*/transformed2);

        // Finally, compute and store the value
        vals.push_back(P_i * P_j);
      }
}



void Dubiner::dp(unsigned d,
                 const mpfr_class & xi,
                 const mpfr_class & eta,
                 std::vector<Point<mpfr_class> > & vals)
{
  // Make sure there are no old values in the vector.
  vals.clear();

  for (unsigned i=0; i<=d; ++i)
    for (unsigned j=0; j<=d; ++j)
      {
        // If degree is too high, continue
        if (i+j > d)
          continue;

        // Map (xi, eta) to barycentric coordinates
        mpfr_class
          zeta0 = xi,
          zeta1 = eta,
          zeta2 = 1 - xi - eta;

        // Compute transformed coordinate for evaluating P_i
        mpfr_class transformed1 = (zeta1-zeta0) / (zeta1+zeta0);

        // Compute P_i^(0,0)
        mpfr_class Pi = this->jacobi(/*n=*/i, /*alpha=*/0, /*beta=*/0, /*x=*/transformed1);

        // Compute d/dx P_i
        mpfr_class dPi_dx = this->djacobi(/*n=*/i, /*alpha=*/0, /*beta=*/0, /*x=*/transformed1);

        // Debugging:
        // std::cout << "d/dx P_{i=" << i << "} = " << dPi_dx << std::endl;

        // Compute d/d(xi) P_i
        mpfr_class dPi_dxi = -2*eta/(xi + eta)/(xi + eta) * dPi_dx;

        // Compute d/d(eta) P_i
        mpfr_class dPi_deta = 2*xi/(xi + eta)/(xi + eta) * dPi_dx;

        // Compute the "scaling" term
        mpfr_class scaling_term = pow(zeta1 + zeta0, i);

        // Compute derivative of scaling term
        mpfr_class dscaling_term = i==0 ? mpfr_class(0.) : mpfr_class(i) * pow(zeta1 + zeta0, i-1);

        // Compute transformed coordinate for evaluating P_j
        mpfr_class transformed2 = 2.*zeta2 - 1.;

        // Compute P_j^(2*i+1,0)
        mpfr_class Pj = this->jacobi(/*n=*/j, /*alpha=*/2*i+1, /*beta=*/0, /*x=*/transformed2);

        // Compute d/d(xi) P_j = d/d(eta) P_j = (-2) * d(P_j)/dx
        mpfr_class dPj = (-2.) * this->djacobi(/*n=*/j, /*alpha=*/2*i+1, /*beta=*/0, /*x=*/transformed2);

        // Debugging (undo dPj scaling):
        // std::cout << "d/dx P_{j=" << j << "}, alpha=" << 2*i+1 << " = " << -0.5*dPj << std::endl;

        // Finally, compute and store the derivatives
        vals.push_back(Point<mpfr_class>(Pi*scaling_term*dPj + Pj*(Pi*dscaling_term + scaling_term*dPi_dxi),
                                         Pi*scaling_term*dPj + Pj*(Pi*dscaling_term + scaling_term*dPi_deta)));
      }
}



void Dubiner::build_H1_projection_matrix(unsigned d,
                                         Matrix<mpfr_class> & matrix)
{
  // The number of polynomials in the Dubiner basis of degree 'd'
  const unsigned Np = (d+1)*(d+2)/2;

  // Resize the input matrix.  Clear it too?
  matrix.resize(Np, Np);

  // We use conical product quadrature of order 2*d to construct the
  // H1 projection matrix.
  Conical conical;
  conical.rule2D(2*d);
  const std::vector<Point<mpfr_class> > & conical_rule_points = conical.get_points();
  const std::vector<mpfr_class> & conical_rule_weights = conical.get_weights();

  // Nq is the number of points in the quadrature rule
  const unsigned Nq = conical_rule_points.size();

  // Fill in the diagonal entries.  We've already independently verified that
  // the off-diagonal entries of the mass matrix are zero.
  for (unsigned q=0; q<Nq; ++q)
    {
      std::vector<mpfr_class> current_vals;

      // Evaluate the Dubiner polynomials at the current qp
      this->p_numeric(d,
                      /*xi=*/  conical_rule_points[q](0),
                      /*eta=*/ conical_rule_points[q](1),
                      current_vals);

      // Evaluate the Dubiner polynomial derivatives at the current qp
      std::vector<Point<mpfr_class> > current_derivs;
      this->dp(d,
               /*xi=*/  conical_rule_points[q](0),
               /*eta=*/ conical_rule_points[q](1),
               current_derivs);

      // Accumulate the mass and Laplace matrix entries.  TODO: only
      // fill in the lower triangle, then copy to the upper triangle.
      for (unsigned i=0; i<Np; ++i)
        for (unsigned j=0; j<Np; ++j)
          {
            matrix(i,j) +=
              conical_rule_weights[q] *
              (current_vals[i]*current_vals[j] +
               current_derivs[i](0)*current_derivs[j](0) +
               current_derivs[i](1)*current_derivs[j](1));
          }
    }
}



mpfr_class Dubiner::jacobi(unsigned n, unsigned alpha, unsigned beta, mpfr_class x)
{
  // We are using the Wikipedia summation formula rather than the
  // 3-term recursion formula. No attempt is made to evaluate
  // factorials carefully, instead we are just using mpq_class objects
  // for this...
  mpfr_class result = 0.;
  for (unsigned s=0; s<=n; ++s)
    {
      // The coefficient of each term is:
      //      (n+alpha)! (n+beta)!
      // --------------------------------
      // (n+alpha-s)! (beta+s)! s! (n-s)!
      mpq_class coeff = 1;
      coeff *= factorial(n+alpha);
      coeff /= factorial(n+alpha-s);
      coeff *= factorial(n+beta);
      coeff /= factorial(beta+s);
      coeff /= factorial(s);
      coeff /= factorial(n-s);

      result += coeff * pow(0.5*(x-1), n-s) * pow(0.5*(x+1), s);
    }

  return result;
}



mpfr_class Dubiner::djacobi(unsigned n, unsigned alpha, unsigned beta, mpfr_class x)
{
  // We are using the Wikipedia summation formula rather than the
  // 3-term recursion formula. No attempt is made to evaluate
  // factorials carefully, instead we are just using mpq_class objects
  // for this...
  mpfr_class result = 0.;
  for (unsigned s=0; s<=n; ++s)
    {
      // The coefficient of each term is:
      //      (n+alpha)! (n+beta)!
      // --------------------------------
      // (n+alpha-s)! (beta+s)! s! (n-s)!
      mpq_class coeff = 1;
      coeff *= factorial(n+alpha);
      coeff /= factorial(n+alpha-s);
      coeff *= factorial(n+beta);
      coeff /= factorial(beta+s);
      coeff /= factorial(s);
      coeff /= factorial(n-s);

      mpfr_class x_term = 0.;
      if (s != n)
        x_term += mpfr_class(0.5)*(n-s) * pow(0.5*(x-1), n-s-1) * pow(0.5*(x+1), s);
      if (s != 0)
        x_term += pow(0.5*(x-1), n-s) * mpfr_class(0.5) * s * pow(0.5*(x+1), s-1);

      result += coeff * x_term;
    }

  return result;
}
