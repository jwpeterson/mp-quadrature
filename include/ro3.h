#ifndef RO3_H
#define RO3_H

// Local includes
#include "matrix.h"
#include "exact.h"
#include "vect.h"

// Support library headers
#include "mpfr.h"
#include "gmpfrxx.h"

// C++ includes
#include <vector>
#include <utility> // std::pair

// A struct which lets the user specify the parameters of an Ro3
// quadrature rule as generally as possible.
struct Ro3
{
  // Constructor. Checks that the parameters are consistent with the dimension.
  Ro3(unsigned int d_in, unsigned int nc_in, unsigned int nv_in,
      unsigned int ne_in, unsigned int ng_in);

  // The 5 "special" functions can be explicitly defaulted for this class.
  Ro3 (const Ro3 &) = default;
  Ro3 (Ro3 &&) = default;
  Ro3 & operator= (const Ro3 &) = default;
  Ro3 & operator= (Ro3 &&) = default;
  ~Ro3() = default;

  // Default ctor: initializes the struct for a d==1 rule with 1 point.
  // Ro3() : d(1), nc(1), nv(0), ne(0), ng(0) {}

  // Enumerate the different types of Orbits.
  enum Orbit : int
    {
      CENTROID,
      VERTEX,
      EDGE,
      GENERAL
    };

  // Returns the index of the first dof in the residual vector
  // corresponding to the type of Orbit, or dim() if there are no
  // orbits of the requested type.
  unsigned int begin(Orbit orb);

  // Returns one past the index of the last dof in the residual vector
  // corresponding to the type of Orbit, or dim() if there are no
  // orbits of the requested type.
  unsigned int end(Orbit orb);

  // Fill the passed-in vectors with upper and lower bounds for the current rule.
  void bounds(std::vector<double> & lb, std::vector<double> & ub);

  // Fill the passed-in vector with a random initial guess that respects the bounds.
  void guess(std::vector<double> & x);

  // Fill the passed-in vector with the indices of the x-dof of each general orbit.
  // This is used by nlopt when setting inequality constraints.
  void inequality_constraint_indices(std::vector<unsigned int> & indices);

  // Compute the residual and Jacobian at u.
  void residual_and_jacobian (std::vector<mpfr_class> * r,
                              Matrix<mpfr_class> * jac,
                              const std::vector<mpfr_class> & u);

  // Returns true if the trial_u solution is "feasible" (i.e. satisfies
  // constraints, false otherwise.
  bool check_feasibility (const std::vector<mpfr_class> & trial_u);

  // The dimension of the space of Ro3-invariant polynomials of degree d.
  unsigned int dim() { return (d*d + 3*d + 6)/6; }

  // The number of quadrature points implied by the current (nc, nv,
  // ne, ng) values.
  unsigned int n_qp() { return nc + 3*(nv + ne + ng); }

  // Check that the current (nc, nv, ne, ng) values are consistent with
  // the current value of d.
  bool check() { return dim() == nc + nv + 2*ne + 3*ng; }

  // Polynomial degree to be integrated exactly
  unsigned int d;

  // Number of centroid orbits (max 1)
  unsigned int nc;

  // Number of vertex orbits (max 1)
  unsigned int nv;

  // Number of edge orbits (max N/2)
  unsigned int ne;

  // Number of "general" orbits (max N/3)
  unsigned int ng;

  // List of monomial exponents for each order that forms a basis for Ro3(d).
  std::vector<std::pair<unsigned int, unsigned int>> polys;

  // A handy multi-precision constant
  mpfr_class one_third;
};

#endif
