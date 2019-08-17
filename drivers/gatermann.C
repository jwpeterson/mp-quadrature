// Functions for computing exact integrals of monomials
// on the reference triangle.
#include "exact.h"
#include "matrix.h"

// C++ includes
#include <vector>
#include <utility>
#include <iostream>

// Function used to (possibly) compute the residual and Jacobian at
// the input u.  If either "r" or "jac" is nullptr, their computation
// is skipped.
void residual_and_jacobian(std::vector<mpfr_class> * r,
                           Matrix<mpfr_class> * jac,
                           const std::vector<mpfr_class> & u);

// A degree 7 rule with 12 points.  This rule can be found in:
//
// K. Gatermann, The construction of symmetric cubature
// formulas for the square and the triangle, Computing 40
// (1988), 229--240.
//
// This rule, which is provably minimal in the number of
// integration points, is said to be 'Ro3 invariant' which
// means that a given set of barycentric coordinates
// (z1,z2,z3) implies the quadrature points (z1,z2),
// (z3,z1), (z2,z3) which are formed by taking the first
// two entries in cyclic permutations of the barycentric
// point.  Barycentric coordinates are related in the
// sense that: z3 = 1 - z1 - z2.
int main()
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // As described in Gatermann, we solve a nonlinear set of equations
  // such that the following set of 12 Ro3-invariant polynomials are
  // integrated exactly:
  // 01.) 1,
  // 02.) x^2,
  // 03.) x^3,
  // 04.) x^2 * y,
  // 05.) x^4,
  // 06.) x^5,
  // 07.) x^4 * y,
  // 08.) x^6,
  // 09.) x^5 * y,
  // 10.) x^7,
  // 11.) x^6 * y,
  // 12.) x^4 * y^2
  std::vector<std::pair<unsigned int, unsigned int>> polys =
    {
      {0,0}, {2,0}, {3,0}, {2,1}, {4,0}, {5,0},
      {4,1}, {6,0}, {5,1}, {7,0}, {6,1}, {4,2}
    };

  for (const auto & p : polys)
    std::cout << "(" << p.first << "," << p.second << ") = "
              << exact_tri(p.first, p.second)
              << std::endl;

  // The nonlinear system of equations is to be solved by Newton's
  // method, so we need to form a residual and Jacobian. The residual is
  // R_i := Q(p_i) - I(p_i)
  // for i=1..12 and p_i one of the reference polynomials from the
  // list above.

  // The quadrature rule consists of four sets of three Ro3-invariant
  // points given in terms of the weight w and the real-valued
  // coordiantes (x, y) of the form:
  // (x, y), (1-x-y, x), (y, 1-x-y)
  // Therefore, the full set of unknowns is:
  // (w1, x1, y1, w2, x2, y2, w3, x3, y3, w4, x4, y4)
  // and the quadrature formula can be written as
  // Q(f) = \sum_{i=1}^4 w_i * (f(xi, yi) + f(1-xi-yi, xi) + f(yi, 1-xi-yi))
  //

  // The solution is initialized to the known solution with 16 digits
  // as given in the paper by Gatermann.
  std::vector<mpfr_class> u =
    {
      2.6517028157436251e-02, // w1
      6.2382265094402118e-02, // x1
      6.7517867073916085e-02, // y1

      4.3881408714446055e-02, // w2
      5.5225456656926611e-02, // x2
      3.2150249385198182e-01, // y2

      2.8775042784981585e-02, // w3
      3.4324302945097146e-02, // x3
      6.6094919618673565e-01, // y3

      6.7493187009802774e-02, // w4
      5.1584233435359177e-01, // x4
      2.7771616697639178e-01, // y4
    };

  std::cout << "u=" << std::endl;
  for (const auto & val : u)
    std::cout << val << std::endl;

  // Compute the residual
  std::vector<mpfr_class> r;
  residual_and_jacobian(&r, nullptr, u);
  // residual_and_jacobian(nullptr, nullptr, u); // test that we can "do nothing"

  // Print the result
  std::cout << "r=" << std::endl;
  for (const auto & val : r)
    std::cout << val << std::endl;

  return 0;
}

void
residual_and_jacobian(std::vector<mpfr_class> * r,
                      Matrix<mpfr_class> * jac,
                      const std::vector<mpfr_class> & u)
{
  // The list of monomial exponents that we are going to test.
  std::vector<std::pair<unsigned int, unsigned int>> polys =
    {
      {0,0}, {2,0}, {3,0}, {2,1}, {4,0}, {5,0},
      {4,1}, {6,0}, {5,1}, {7,0}, {6,1}, {4,2}
    };

  // Allocate space for residual vector
  if (r)
    r->resize(polys.size());

  // Loop over (w,x,y) triples in u
  for (unsigned int i=0; i<polys.size(); i++)
    {
      // The powers of x and y in the monomial we are currently integrating.
      unsigned int xpower = polys[i].first;
      unsigned int ypower = polys[i].second;

      // Accumlate quadrature rule contribution for r(i).
      for (unsigned int q=0; q<u.size(); q+=3)
        {
          // The unknowns are ordered in terms of (w,x,y) triples.
          mpfr_class w = u[q];
          mpfr_class x = u[q+1];
          mpfr_class y = u[q+2];

          // The implied third barycentric coordinate
          mpfr_class z = mpfr_class(1) - x - y;

          if (r)
            (*r)[i] += w * (pow(x, xpower) * pow(y, ypower) +
                            pow(z, xpower) * pow(x, ypower) +
                            pow(y, xpower) * pow(z, ypower));
        }

      // Subtract off the true integral value, I(p_i)
      if (r)
        (*r)[i] -= exact_tri(xpower, ypower);
    }
}
