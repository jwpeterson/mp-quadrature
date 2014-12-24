#include <cstdlib> // std::abort
#include "dubiner.h"
#include "common_definitions.h"

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




void Dubiner::p(unsigned d,
                const mpfr_class & xi,
                const mpfr_class & eta,
                std::vector<mpfr_class> & vals)
{
  // Make sure there are no old values in the vector.
  vals.clear();

  // Compute the barycenter coordinate values
  mpfr_class
    zeta0 = xi,
    zeta1 = eta,
    zeta2 = 1. - xi - eta;

  /* (0,0) */ vals.push_back(1);

  if (d >= 1)
    {
      /* (0,1) */ vals.push_back(3*zeta2 - 1);
      /* (1,0) */ vals.push_back(-zeta0 + zeta1);
    }
  if (d >= 2)
    {
      // Note: gmpfrxx does not like it when you have an "L" suffix on your constants!
      // We need to be careful here if there are fractions that do not have a finite
      // decimal representation -- they will not be full precision unless you explicitly
      // cast stuff to the mpfr_class type...
      /* (0,2) */ vals.push_back(3*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (3.0/2.0)*pow(2*zeta2 - 1, 2) - 3.0/2.0);
      /* (1,1) */ vals.push_back((-zeta0 + zeta1)*(5*zeta2 - 1));
      /* (2,0) */ vals.push_back(pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2));
    }
  if (d >= 3)
    {
      /* (0,3) */ vals.push_back(4*pow(zeta2, 3) + 9*pow(zeta2, 2)*(2*zeta2 - 2) + 3*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3));
      /* (1,2) */ vals.push_back((-zeta0 + zeta1)*(10*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (5.0/2.0)*pow(2*zeta2 - 1, 2) - 5.0/2.0));
      /* (2,1) */ vals.push_back((7*zeta2 - 1)*(pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2)));
      /* (3,0) */ vals.push_back(-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3));
    }
  if (d >= 4)
    {
      /* (0,4) */ vals.push_back(5*pow(zeta2, 4) + 20*pow(zeta2, 3)*(2*zeta2 - 2) + 15*pow(zeta2, 2)*pow(2*zeta2 - 2, 2) + (5.0/2.0)*zeta2*pow(2*zeta2 - 2, 3) + (1.0/16.0)*pow(2*zeta2 - 2, 4));
      /* (1,3) */ vals.push_back((-zeta0 + zeta1)*(20*pow(zeta2, 3) + (45.0/2.0)*pow(zeta2, 2)*(2*zeta2 - 2) + (9.0/2.0)*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3)));
      /* (2,2) */ vals.push_back((pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2))*(21*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (7.0/2.0)*pow(2*zeta2 - 1, 2) - 7.0/2.0));
      /* (3,1) */ vals.push_back((9*zeta2 - 1)*(-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3)));
      /* (4,0) */ vals.push_back(pow(zeta0, 4) - 16*pow(zeta0, 3)*zeta1 + 36*pow(zeta0, 2)*pow(zeta1, 2) - 16*zeta0*pow(zeta1, 3) + pow(zeta1, 4));
    }
  if (d >= 5)
    this->dubiner_5th(zeta0, zeta1, zeta2, vals);
  if (d >= 6)
    this->dubiner_6th(zeta0, zeta1, zeta2, vals);
  if (d >= 7)
    this->dubiner_7th(zeta0, zeta1, zeta2, vals);
  if (d >= 8)
    this->dubiner_8th(zeta0, zeta1, zeta2, vals);
  if (d >= 9)
    this->dubiner_9th(zeta0, zeta1, zeta2, vals);
  if (d >= 10)
    this->dubiner_10th(zeta0, zeta1, zeta2, vals);
  if (d >= 11)
    this->dubiner_11th(zeta0, zeta1, zeta2, vals);
  if (d >= 12)
    this->dubiner_12th(zeta0, zeta1, zeta2, vals);
  if (d >= 13)
    this->dubiner_13th(zeta0, zeta1, zeta2, vals);
  if (d >= 14)
    this->dubiner_14th(zeta0, zeta1, zeta2, vals);
  if (d >= 15)
    this->dubiner_15th(zeta0, zeta1, zeta2, vals);
  if (d >= 16)
    this->dubiner_16th(zeta0, zeta1, zeta2, vals);
  if (d >= 17)
    this->dubiner_17th(zeta0, zeta1, zeta2, vals);
  if (d >= 18)
    this->dubiner_18th(zeta0, zeta1, zeta2, vals);
  if (d >= 19)
    this->dubiner_19th(zeta0, zeta1, zeta2, vals);
  if (d >= 20)
    this->dubiner_20th(zeta0, zeta1, zeta2, vals);
  if (d >= 21)
    this->dubiner_21st(zeta0, zeta1, zeta2, vals);
  if (d >= 22)
    this->dubiner_22nd(zeta0, zeta1, zeta2, vals);
  if (d >= 23)
    this->dubiner_23rd(zeta0, zeta1, zeta2, vals);
  if (d >= 24)
    this->dubiner_24th(zeta0, zeta1, zeta2, vals);
  if (d >= 25)
    this->dubiner_25th(zeta0, zeta1, zeta2, vals);
  if (d >= 26)
    this->dubiner_26th(zeta0, zeta1, zeta2, vals);
  if (d >= 27)
    this->dubiner_27th(zeta0, zeta1, zeta2, vals);
  if (d >= 28)
    this->dubiner_28th(zeta0, zeta1, zeta2, vals);
  if (d >= 29)
    this->dubiner_29th(zeta0, zeta1, zeta2, vals);
  if (d >= 30)
    this->dubiner_30th(zeta0, zeta1, zeta2, vals);
}



void Dubiner::build_H1_projection_matrix(unsigned d,
                                         Matrix<mpfr_class> & matrix)
{
  // Currently only works for d==10
  if (d != 10)
    {
      std::cerr << "Can only build H1-projection matrix for d==10!" << std::endl;
      std::abort();
    }

  // Get the "orthogonality coeffs", which are int(phi_i * phi_j) for
  // the requested value of d.
  std::vector<mpq_class> coeffs;
  this->orthogonality_coeffs(d, coeffs);

  // The number of rows/cols in the matrix
  const unsigned N = coeffs.size();

  // Resize the input matrix.
  matrix.resize(coeffs.size(), coeffs.size());

  // Fill in the diagonal entries
  for (unsigned i=0; i<N; ++i)
    matrix(i,i) = coeffs[i];

  // Fill in the off-diagonal parts of the matrix.  These were
  // pre-computed in Python...
  for (unsigned i=0; i<N; ++i)
    for (unsigned j=0; j<N; ++j)
      {
        // The upper triangle is empty
        if (j > i)
          matrix(i,j) += laplace_matrix[j][i];
        else
          matrix(i,j) += laplace_matrix[i][j];
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

// Local Variables:
// truncate-lines: t
// End:
