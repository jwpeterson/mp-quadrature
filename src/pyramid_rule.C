#include "pyramid_rule.h"

void PyramidRule::generate(unsigned int order)
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

  // Compute the Gauss rule points/weights
  std::vector<mpfr_class> gauss_x, gauss_w;
  gauss_rule(n_points, gauss_x, gauss_w);

  // Compute a Jacobi rule with alpha==2.  The Jacobi weights by
  // default sum to 1 and lie in [-1,1], we want the weights to sum to
  // 1/3 and lie in [0,1] for this application.
  Jacobi p1(/*alpha=*/2., /*beta=*/0.);
  p1.rule(n_points);
  p1.scale_weights(mpfr_class(1)/3);
  p1.scale_points(0., 1.);

  // Get const references to the jacobi points and weights vectors.
  const std::vector<mpfr_class> & jacobi_x = p1.get_points();
  const std::vector<mpfr_class> & jacobi_w = p1.get_weights();

  // Compute the conical product rule, with space for n^3 entries.
  // See also the code in LibMesh's src/quadrature/quadrature_conical.C
  x.resize(n_points*n_points*n_points);
  w.resize(n_points*n_points*n_points);
  unsigned int gp = 0;
  for (unsigned int i=0; i<n_points; i++)
    for (unsigned int j=0; j<n_points; j++)
      for (unsigned int k=0; k<n_points; k++)
        {
          // Note: Access the 1D arrays from [1] ... [n]
          mpfr_class zk = jacobi_x[k+1];

          x[gp](0) = (1-zk) * gauss_x[i+1];
          x[gp](1) = (1-zk) * gauss_x[j+1];
          x[gp](2) = zk;
          w[gp] = gauss_w[i+1] * gauss_w[j+1] * jacobi_w[k+1];
          gp++;
        }
}
