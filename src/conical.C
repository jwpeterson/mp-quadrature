#include "conical.h"

void Conical::rule2D(unsigned int order)
{
  // Clear any existing data
  x.clear();
  w.clear();

  // If the user asks for an even order, he gets the next highest odd order for free
  if (order % 2 == 0)
    order++;

  // The number of points in each of the constituent 1D rules.
  // In 1D, order = 2*n-1, so n = (order+1)/2
  unsigned n_points = (order+1)/2;

  // Compute the 'n_points' Gauss rule points/weights
  std::vector<mpfr_class> gauss_x, gauss_w;
  gauss_rule(n_points, gauss_x, gauss_w);

  // Scale the Gauss points so they lie in [0,1].  We should probably make
  // a Gauss class and have this be a member. See eg: Jacobi
  mpfr_class zero(0.0), one(1.0);
  {
    mpfr_class a (0.5*(one-zero));
    mpfr_class b (0.5*(zero+one));
    for (unsigned int j=1; j<gauss_x.size(); ++j)
      {
        gauss_x[j] = a*gauss_x[j] + b;
        gauss_w[j] *= 0.5;
      }
  }

  // Compute the 'n_points' Jacobi rule points/weights for alpha=1 and
  // scale the points and weights.
  Jacobi jacobi(/*alpha=*/1., /*beta=*/0.);
  jacobi.rule(n_points);
  jacobi.scale_weights(0.5);
  jacobi.scale_points(zero, one);

  // Get const references to the jacobi points and weights vectors.
  const std::vector<mpfr_class>& jacobi_x = jacobi.get_points();
  const std::vector<mpfr_class>& jacobi_w = jacobi.get_weights();

  // Compute the conical product rule with n^2 entries.
  x.resize(n_points*n_points);
  w.resize(n_points*n_points);
  unsigned int gp = 0;
  for (unsigned int i=0; i<n_points; i++)
    for (unsigned int j=0; j<n_points; j++)
      {
        // Note: Access the Jacobi and Gauss 1D arrays from [1] ... [n]
        x[gp](0) = jacobi_x[j+1];                     //s[j];
        x[gp](1) = gauss_x[i+1] * (1.-jacobi_x[j+1]); //r[i]*(1.-s[j]);
        w[gp]    = gauss_w[i+1] * jacobi_w[j+1];      //A[i]*B[j];
        gp++;
      }

}


