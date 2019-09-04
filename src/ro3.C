#include "ro3.h"


Ro3::Ro3(unsigned int d_in, unsigned int nc_in, unsigned int nv_in,
         unsigned int ne_in, unsigned int ng_in) :
  d(d_in), nc(nc_in), nv(nv_in), ne(ne_in), ng(ng_in)
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
    {
      {0,0},                                         // const
      {2,0},                                         // 2nd
      {3,0},  {2,1},                                 // 3rd
      {4,0},                                         // 4th
      {5,0},  {4,1},                                 // 5th
      {6,0},  {5,1},  {4,2},                         // 6th
      {7,0},  {6,1},                                 // 7th
      {8,0},  {7,1},  {6,2},                         // 8th
      {9,0},  {8,1},  {7,2},  {6,3},                 // 9th
      {10,0}, {9,1},  {8,2},                         // 10th
      {11,0}, {10,1}, {9,2},  {8,3},                 // 11th
      {12,0}, {11,1}, {10,2}, {9,3},  {8,4},         // 12th
      {13,0}, {12,1}, {11,2}, {10,3},                // 13th
      {14,0}, {13,1}, {12,2}, {11,3}, {10,4},        // 14th
      {15,0}, {14,1}, {13,2}, {12,3}, {11,4}, {10,5} // 15th
    };

  // Throw an error if we have not tabulated enough basis polynomials yet.
  if (polys.size() < dim())
    {
      std::cout << "Not enough polynomials for d = " << d << std::endl;
      exit(1);
    }
}



unsigned int Ro3::first_dof(Orbit orb)
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



unsigned int Ro3::last_dof(Orbit orb)
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
  for (unsigned int i=first_dof(CENTROID); i<last_dof(CENTROID); ++i)
    ub[i] = 0.5;

  // Vertex orbits have a weight dof that appears three times.
  for (unsigned int i=first_dof(VERTEX); i<last_dof(VERTEX); ++i)
    ub[i] = 1. / 6;

  // Edge orbits have a weight dof (appearing 3 times) and a spatial dof.
  for (unsigned int i=first_dof(EDGE); i<last_dof(EDGE); i+=2)
    {
      ub[i] = 1. / 6;
      ub[i+1] = 1.;
    }

  // General orbits have a weight dof and two spatial dofs.
  for (unsigned int i=first_dof(GENERAL); i<last_dof(GENERAL); i+=3)
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
  for (unsigned int i=first_dof(CENTROID); i<last_dof(CENTROID); ++i)
    x[i] = 0.5 * double(random())/RAND_MAX;

  // Vertex orbits have a weight dof that appears three times.
  for (unsigned int i=first_dof(VERTEX); i<last_dof(VERTEX); ++i)
    x[i] = 1. / 6 * double(random())/RAND_MAX;

  // Edge orbits have a weight dof (appearing 3 times) and a spatial dof.
  for (unsigned int i=first_dof(EDGE); i<last_dof(EDGE); i+=2)
    {
      x[i] = 1. / 6 * double(random())/RAND_MAX;
      x[i+1] = 1. * double(random())/RAND_MAX;
    }

  // General orbits have a weight dof and two spatial dofs.
  for (unsigned int i=first_dof(GENERAL); i<last_dof(GENERAL); i+=3)
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
  for (unsigned int i=first_dof(GENERAL); i<last_dof(GENERAL); i+=3)
    indices.push_back(i+1);
}
