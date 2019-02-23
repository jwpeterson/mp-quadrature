#include <cstdlib> // std::abort
#include "dubiner.h"
#include "common_definitions.h"
#include "conical.h"

void Dubiner::p(unsigned d,
                const mpfr_class & xi,
                const mpfr_class & eta,
                std::vector<mpfr_class> & vals,
                std::vector<Point<mpfr_class> > & gradients)
{
  // Make sure there are no old values in the vectors.
  vals.clear();
  gradients.clear();

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

        // Compute P_i^(0,0) and derivative
        std::pair<mpfr_class, mpfr_class> jacobi_i = this->jacobi(/*n=*/i, /*alpha=*/0, /*beta=*/0, /*x=*/transformed1);

        // Shortcut names
        const mpfr_class & P_i = jacobi_i.first;
        const mpfr_class & dPi_dx = jacobi_i.second;

        // Compute the "scaling" term
        mpfr_class scaling_term = pow(zeta1 + zeta0, i);

        // Compute transformed coordinate for evaluating P_j
        mpfr_class transformed2 = 2.*zeta2 - 1.;

        // Compute P_j^(2*i+1,0) and derivative
        std::pair<mpfr_class, mpfr_class> jacobi_j = this->jacobi(/*n=*/j, /*alpha=*/2*i+1, /*beta=*/0, /*x=*/transformed2);

        // Shortcut names
        const mpfr_class & P_j = jacobi_j.first;
        const mpfr_class & dPj_dx = jacobi_j.second;

        // Compute and store the value
        vals.push_back(P_i * scaling_term * P_j);

        // Compute derivatives

        // Compute d/d(xi) P_i
        mpfr_class dPi_dxi = -2*eta/(xi + eta)/(xi + eta) * dPi_dx;

        // Compute d/d(eta) P_i
        mpfr_class dPi_deta = 2*xi/(xi + eta)/(xi + eta) * dPi_dx;

        // Compute derivative of scaling term
        mpfr_class dscaling_term = i==0 ? mpfr_class(0.) : mpfr_class(i) * pow(zeta1 + zeta0, i-1);

        // Compute d/d(xi) P_j = d/d(eta) P_j = (-2) * d(P_j)/dx
        mpfr_class dPj = (-2.) * dPj_dx;

        // Finally, compute and store the derivatives
        gradients.push_back(Point<mpfr_class>(P_i*scaling_term*dPj + P_j*(P_i*dscaling_term + scaling_term*dPi_dxi),
                                              P_i*scaling_term*dPj + P_j*(P_i*dscaling_term + scaling_term*dPi_deta)));
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
      std::vector<Point<mpfr_class> > current_derivs;

      // Evaluate the Dubiner polynomials at the current qp
      this->p(d,
              /*xi=*/  conical_rule_points[q](0),
              /*eta=*/ conical_rule_points[q](1),
              current_vals,
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



std::pair<mpfr_class, mpfr_class> Dubiner::jacobi(unsigned n, unsigned alpha, unsigned beta, mpfr_class x)
{
  // We are using the Wikipedia summation formula rather than the
  // 3-term recursion formula. No attempt is made to evaluate
  // factorials carefully, instead we are just using mpq_class objects
  // for this...
  mpfr_class
    val = 0.,
    deriv = 0.;

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

      // Compute the value
      val += coeff * pow(0.5*(x-1), n-s) * pow(0.5*(x+1), s);

      // Compute the derivative
      mpfr_class x_term = 0.;
      if (s != n)
        x_term += mpfr_class(0.5)*(n-s) * pow(0.5*(x-1), n-s-1) * pow(0.5*(x+1), s);
      if (s != 0)
        x_term += pow(0.5*(x-1), n-s) * mpfr_class(0.5) * s * pow(0.5*(x+1), s-1);

      deriv += coeff * x_term;
    }

  return std::make_pair(val, deriv);
}

std::pair<double, double>
Dubiner::jacobi(unsigned n, unsigned alpha, unsigned beta, double x)
{
  if (n == 0)
    return std::make_pair(1,0);

  Real p0 = 1;
  Real p1 = (alpha + 1) + (alpha + beta + 2) * 0.5 * (x - 1);

  // The recurrence relation we're using only requires one previous
  // derivative value, but it is from two iterations ago, so we still
  // need to keep track of two values.
  Real dp0 = 0;
  Real dp1 = 0.5 * (alpha + beta + 2);

  unsigned int i = 1;
  while (i < n)
    {
      // Swapping saves a temporary, since this immediately updates p0
      // and then p1 is updated on the next line. Note that p0 and p1
      // appear in opposite positions than is usual in the update
      // formula because of this swap!
      std::swap(p0, p1);
      p1 = ((2*n + alpha + beta - 1) *
        ((2*n + alpha + beta) * (2*n + alpha + beta - 2) * x + alpha * alpha - beta * beta) * p0
        - 2 * (n + alpha - 1) * (n + beta - 1) * (2*n + alpha + beta) * p1) /
        (2*n * (n + alpha + beta) * (2*n + alpha + beta - 2));

      // Note that dp1 appears in the formula below, but because we
      // swapped, it's really dp0. Also, we have already updated the
      // values, so:
      // p1 is now P_{i+1}(x)
      // p0 is now P_{i}(x)
      std::swap(dp0, dp1);
      dp1 = (2*i + 1) * p0 + dp1;
      ++i;
    }
  return std::make_pair(p1, dp1);
}

void
Dubiner::compare_jacobi()
{
  unsigned int alpha = 1;
  unsigned int beta = 0;
  unsigned int n = 2;

  double x = 0.2;
  mpfr_class mp_x = mpfr_class(2.) / mpfr_class(10.);
  auto mp_result = this->jacobi(n, alpha, beta, mp_x);
  auto double_result = this->jacobi(n, alpha, beta, x);
  std::cout << "mp result: " << mp_result.first
            << ", double result: " << double_result.first
            << std::endl;
}
