// Functions for computing exact integrals of monomials
// on the reference triangle.
#include "exact.h"

// C++ includes
#include <vector>
#include <utility>
#include <iostream>

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

  return 0;
}
