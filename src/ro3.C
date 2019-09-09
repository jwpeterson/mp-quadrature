#include "ro3.h"


Ro3::Ro3(unsigned int d_in, unsigned int nc_in, unsigned int nv_in,
         unsigned int ne_in, unsigned int ng_in) :
  d(d_in),
  nc(nc_in),
  nv(nv_in),
  ne(ne_in),
  ng(ng_in),
  one_third(mpfr_class(1) / 3)
{
  if (nc > 1 || nv > 1 || !check())
    {
      std::cerr << "Parameters ("
                << nc << ","
                << nv << ","
                << ne << ","
                << ng << ")"
                << " are inconsistent with degree " << d
                << std::endl;
      exit(1);
    }

  // Initialize the polynomial exponents array.
  polys =
    {                                                         //        a_d
      {0,0},                                                  // const  1
      {2,0},                                                  // 2nd    1
      {3,0},  {2,1},                                          // 3rd    2]
      {4,0},                                                  // 4th    1]
      {5,0},  {4,1},                                          // 5th    2]
      {6,0},  {5,1},  {4,2},                                  // 6th    3
      {7,0},  {6,1},                                          // 7th    2
      {8,0},  {7,1},  {6,2},                                  // 8th    3
      {9,0},  {8,1},  {7,2},  {6,3},                          // 9th    4]
      {10,0}, {9,1},  {8,2},                                  // 10th   3]
      {11,0}, {10,1}, {9,2},  {8,3},                          // 11th   4]
      {12,0}, {11,1}, {10,2}, {9,3},  {8,4},                  // 12th   5
      {13,0}, {12,1}, {11,2}, {10,3},                         // 13th   4
      {14,0}, {13,1}, {12,2}, {11,3}, {10,4},                 // 14th   5
      {15,0}, {14,1}, {13,2}, {12,3}, {11,4}, {10,5},         // 15th   6]
      {16,0}, {15,1}, {14,2}, {13,3}, {12,4},                 // 16th   5]
      {17,0}, {16,1}, {15,2}, {14,3}, {13,4}, {12,5},         // 17th   6]
      {18,0}, {17,1}, {16,2}, {15,3}, {14,4}, {13,5}, {12,6}, // 18th   7
      {19,0}, {18,1}, {17,2}, {16,3}, {15,4}, {14,5},         // 19th   6
      {20,0}, {19,1}, {18,2}, {17,3}, {16,4}, {15,5}, {14,6}  // 20th   7
    };

  // Throw an error if we have not tabulated enough basis polynomials yet.
  if (polys.size() < dim())
    {
      std::cout << "Not enough polynomials for d = " << d << std::endl;
      exit(1);
    }
}



unsigned int Ro3::begin(Orbit orb)
{
  switch (orb)
    {
    case CENTROID:
      return nc ? 0 : dim();
    case VERTEX:
      if (!nv) return dim();
      return nc;
    case EDGE:
      if (!ne) return dim();
      return nc + nv;
    case GENERAL:
      if (!ng) return dim();
      return nc + nv + 2*ne;
    default:
      std::cerr << "Unknown Orbit type." << std::endl;
      exit(1);
    }
}



unsigned int Ro3::end(Orbit orb)
{
  switch (orb)
    {
    case CENTROID:
      return nc ? 1 : dim();
    case VERTEX:
      if (!nv) return dim();
      return nc + nv;
    case EDGE:
      if (!ne) return dim();
      return nc + nv + 2*ne;
    case GENERAL:
      return dim();
    default:
      std::cerr << "Unknown Orbit type." << std::endl;
      exit(1);
    }
}

void Ro3::bounds(std::vector<double> & lb, std::vector<double> & ub)
{
  // Lower bounds are currently all zero.
  lb.clear();
  lb.resize(dim());

  ub.clear();
  ub.resize(dim());

  // Centroid orbits have a single weight dof, and it just needs to
  // be less than the reference element volume.
  for (unsigned int i=begin(CENTROID); i<end(CENTROID); ++i)
    ub[i] = 0.5;

  // Vertex orbits have a weight dof that appears three times.
  for (unsigned int i=begin(VERTEX); i<end(VERTEX); ++i)
    ub[i] = 1. / 6;

  // Edge orbits have a weight dof (appearing 3 times) and a spatial dof.
  for (unsigned int i=begin(EDGE); i<end(EDGE); i+=2)
    {
      ub[i] = 1. / 6;
      ub[i+1] = 1.;
    }

  // General orbits have a weight dof and two spatial dofs.
  for (unsigned int i=begin(GENERAL); i<end(GENERAL); i+=3)
    {
      ub[i] = 1. / 6;
      ub[i+1] = 1.;
      ub[i+2] = 1.;
    }
}

void Ro3::guess(std::vector<double> & x)
{
  x.clear();
  x.resize(dim());

  // Centroid orbits have a single weight dof, and it just needs to
  // be less than the reference element volume.
  for (unsigned int i=begin(CENTROID); i<end(CENTROID); ++i)
    x[i] = 0.5 * double(random())/RAND_MAX;

  // Vertex orbits have a weight dof that appears three times.
  for (unsigned int i=begin(VERTEX); i<end(VERTEX); ++i)
    x[i] = 1. / 6 * double(random())/RAND_MAX;

  // Edge orbits have a weight dof (appearing 3 times) and a spatial dof.
  for (unsigned int i=begin(EDGE); i<end(EDGE); i+=2)
    {
      x[i] = 1. / 6 * double(random())/RAND_MAX;
      x[i+1] = 1. * double(random())/RAND_MAX;
    }

  // General orbits have a weight dof and two spatial dofs.
  for (unsigned int i=begin(GENERAL); i<end(GENERAL); i+=3)
    {
      x[i] = 1. / 6 * double(random())/RAND_MAX;
      x[i+1] = 1. * double(random())/RAND_MAX;
      x[i+2] = 1. * double(random())/RAND_MAX;
    }
}


void Ro3::inequality_constraint_indices(std::vector<unsigned int> & indices)
{
  indices.clear();

  // Push back the index of each general orbit's x-coordinate.
  for (unsigned int i=begin(GENERAL); i<end(GENERAL); i+=3)
    indices.push_back(i+1);
}



void Ro3::residual_and_jacobian (std::vector<mpfr_class> * r,
                                 Matrix<mpfr_class> * jac,
                                 const std::vector<mpfr_class> & u)
{
  // If there's nothing to do, then there's nothing to do.
  if (!r && !jac)
    return;

  // Avoid calling dim() multiple times.
  const unsigned int N = dim();

  // Zero any previous values and allocate space for residual and Jacobian, if required.
  if (r)
    {
      r->clear();
      r->resize(N);
    }
  if (jac)
    {
      jac->clear();
      jac->resize(N, N);
    }

  // Compute residual and Jacobian contributions For each basis function i.
  for (unsigned int i=0; i<N; i++)
    {
      // The powers of x and y in the monomial we are currently integrating.
      unsigned int xpower = polys[i].first;
      unsigned int ypower = polys[i].second;

      // Residual & Jacobian contributions due to centroid point, if any.
      // The evaluation point is fixed at (x,y) = (1/3, 1/3).
      for (unsigned int q=begin(CENTROID); q<end(CENTROID); ++q)
        {
          if (r)
            (*r)[i] += u[q] * pow(one_third, xpower) * pow(one_third, ypower);

          if (jac)
            (*jac)(i,q) += pow(one_third, xpower) * pow(one_third, ypower);
        }

      // Residual & Jacobian contributions due to vertex points, if any.
      for (unsigned int q=begin(VERTEX); q<end(VERTEX); ++q)
        {
          mpfr_class spatial(0);

          if (xpower == 0 && ypower == 0)
            spatial = mpfr_class(3);

          else if ((xpower == 0 && ypower > 0) ||
                   (xpower > 0 && ypower == 0))
            spatial = mpfr_class(1);

          else // both xpower > 0 && ypower > 0
            spatial = 0;

          if (r)
            (*r)[i] += u[q] * spatial;

          if (jac)
            (*jac)(i,q) += spatial;
        }

      // Residual & Jacobian contributions due to edge orbits.
      for (unsigned int q=begin(EDGE); q<end(EDGE); q+=2)
        {
          mpfr_class w = u[q];
          mpfr_class x = u[q+1];
          // mpfr_class y(0);

          // The implied third barycentric coordinate.
          // The three points of this orbit are: (x,0), (z,x), (0,z)
          mpfr_class z = mpfr_class(1) - x; // - y

          // The "spatial" part is needed by both the residual and Jacobian,
          // so we can always compute it.

          // Note: we go to some extra trouble here to absolutely
          // avoid any possible issues with 0^0, even though I believe
          // it should be guaranteed that pow(base, 0) returns 1 for
          // any base. This way there is no possible confusion for the
          // case of 0^0, which we also want to return 1 for.
          mpfr_class spatial(0);

          if (xpower == 0 && ypower == 0)
            spatial = mpfr_class(3);

          else if (xpower == 0 && ypower > 0)
            spatial = pow(x, ypower) + pow(z, ypower);

          else if (xpower > 0 && ypower == 0)
            spatial = pow(x, xpower) + pow(z, xpower);

          else // both xpower > 0 && ypower > 0
            spatial = pow(z, xpower) * pow(x, ypower);

          // Compute residual contribution, if required.
          if (r)
            (*r)[i] += w * spatial;

          // Compute Jacobian contribution, if required.
          if (jac)
            {
              // Derivative wrt w
              (*jac)(i, q) += spatial;

              // Derivative wrt x is zero for the constant polynomial case.
              if (xpower == 0 && ypower == 0)
                (*jac)(i, q+1) += 0;

              else if (xpower == 0 && ypower > 0)
                (*jac)(i, q+1) += w * ypower * (pow(x, ypower-1) - pow(z, ypower-1));

              else if (xpower > 0 && ypower == 0)
                (*jac)(i, q+1) += w * xpower * (pow(x, xpower-1) - pow(z, xpower-1));

              else // both xpower > 0 && ypower > 0
                {
                  unsigned int xpm1 = xpower - 1;
                  unsigned int ypm1 = ypower - 1;

                  // Derivative wrt x. In this case, both
                  // y**ypower = 0
                  // y**xpower = 0
                  // so we have dropped those terms.
                  (*jac)(i, q+1) += w *
                    (mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) +
                     pow(z, xpower) * ypower * pow(x, ypm1));
                }
            }
        }

      // Residual & Jacobian contributions due to general orbits.
      for (unsigned int q=begin(GENERAL); q<end(GENERAL); q+=3)
        {
          // The unknowns are ordered in terms of (w,x,y) triples.
          mpfr_class w = u[q];
          mpfr_class x = u[q+1];
          mpfr_class y = u[q+2];

          // The implied third barycentric coordinate
          mpfr_class z = mpfr_class(1) - x - y;

          // The "spatial" part is needed by both the residual and Jacobian,
          // so we can always compute it.
          mpfr_class spatial =
            pow(x, xpower) * pow(y, ypower) +
            pow(z, xpower) * pow(x, ypower) +
            pow(y, xpower) * pow(z, ypower);

          // Compute residual contribution, if required.
          if (r)
            (*r)[i] += w * spatial;

          // Compute Jacobian contribution, if required.
          if (jac)
            {
              // Derivative wrt w
              (*jac)(i, q) += spatial;

              // We are differentiating polynomials, so if xpower or
              // ypower is 0, the derivative of the corresponding term
              // will be zero. In that case we don't want to subtract
              // 1 from the unsigned variable which represents the
              // exponent, since that would cause it to wrap, possibly
              // leading to other issues. Therefore we will just
              // define the power to still be zero in that case...
              // Note that the -1 terms come from differentiating "z"
              // wrt either x or y.
              unsigned int xpm1 = xpower > 0 ? xpower-1 : 0;
              unsigned int ypm1 = ypower > 0 ? ypower-1 : 0;

              // Derivative wrt x
              (*jac)(i,q+1) += w *
                (xpower * pow(x, xpm1) * pow(y, ypower) +
                 (mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) +
                  pow(z, xpower) * ypower * pow(x, ypm1)) +
                 pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1));

              // Derivative wrt y
              (*jac)(i,q+2) += w *
                (pow(x, xpower) * ypower * pow(y, ypm1) +
                 mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) +
                 (xpower * pow(y, xpm1) * pow(z, ypower) +
                  pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1)));
            }
        }

      // Subtract off the true integral value, I(p_i)
      if (r)
        (*r)[i] -= exact_tri(xpower, ypower);
    } // end loop over i
}



bool Ro3::check_feasibility (const std::vector<mpfr_class> & trial_u)
{
  const unsigned int N = dim();

  if (N != trial_u.size())
    {
      std::cout << "Error, wrong size solution in check_feasibility()" << std::endl;
      exit(1);
    }

  // 1.) No parameters can be negative.
  for (unsigned int q=0; q<trial_u.size(); ++q)
    if (trial_u[q] < 0.)
      {
        // Debugging:
        // std::cout << "trial_u is infeasible due to parameter <= 0!"
        //           << std::endl;
        // print(trial_u);

        return false;
      }

  // 2.) Check if 1-x-y < 0 in general orbits
  for (unsigned int q=begin(GENERAL); q<end(GENERAL); q+=3)
    {
      mpfr_class x = trial_u[q+1];
      mpfr_class y = trial_u[q+2];
      if (mpfr_class(1) - x - y < 0)
        {
          // Debugging:
          // std::cout << "trial_u is infeasible due to 1-x-y <= 0!"
          //           << std::endl;
          // print(trial_u);

          return false;
        }
    }

  // If we made it here without returning, must be feasible!
  return true;
}
