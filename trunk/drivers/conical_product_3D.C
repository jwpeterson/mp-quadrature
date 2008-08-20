#include <iomanip>
#include "gauss.h"
#include "jacobi.h"

// In this program we compute multi-precision
// points and weights to be used in 3D Conical Product
// formulae for tetrahedra.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256); 

  // Read number of points in rule from command line
  unsigned int n=6;
  if (argc > 1)
    n=atoi(argv[1]);
  
  const unsigned int n_total = n*n*n;
  std::cout << "\nComputing 3D conical product rule with n^2=" << n_total << " points." << std::endl;

  // Compute the Gauss rule points/weights
  std::vector<mpfr_class> gauss_x, gauss_w;
  gauss_rule(n, gauss_x, gauss_w);

  // Scale the Gauss points so they lie in [0,1].  We should probably make
  // a Gauss class and have this be a member. See eg: Jacobi
  mpfr_class zero(0.0), one(1.0);
  {
    mpfr_class a ( 0.5*(one-zero) );
    mpfr_class b ( 0.5*(zero+one) );
    for (unsigned int j=1; j<gauss_x.size(); ++j)
      {
	gauss_x[j] = a*gauss_x[j] + b;
	gauss_w[j] *= 0.5;
      }
  }

  // Compute the Jacobi rules' points/weights:
  //
  // p1 := alpha=1, beta=0
  const Real alpha1=1.0, beta1=0.0;
  Jacobi p1(alpha1,beta1);
  p1.rule(n);
  p1.scale_weights(0.5);
  p1.scale_points(zero, one);


  // p2 := alpha=2, beta=0
  const Real alpha2=2.0, beta2=0.0;
  Jacobi p2(alpha2,beta2);
  p2.rule(n);
  p2.scale_weights(0.5);
  p2.scale_points(zero, one);

  
  // Scale Jacobi weights so they sum to 1/3 (alpha==2) or 1/2 (alpha==1)
  if (alpha==2.0)
    {
      mpfr_class one_third(1.0);
      one_third /= 3.0;
      p1.scale_weights(one_third);
    }
  else if (alpha==1.0)
    {
      p1.scale_weights(0.5);
    }
  else
    {
      std::cout << "Warning: weights unscaled!" << std::endl;
    }


  // Print the result
  // p1.printxw();

  // Get const references to the jacobi points and weights vectors.
  const std::vector<mpfr_class>& jacobi1_x = p1.get_points();
  const std::vector<mpfr_class>& jacobi1_w = p1.get_weights();

  // Compute the conical product rule, with space for n^3 entries.
  // See also the code in LibMesh's src/quadrature/quadrature.C
  std::vector<mpfr_class> conical_x(n_total), conical_y(n_total), conical_z(n_total), conical_w(n_total);
  unsigned int gp = 0;
  for (unsigned int i=0; i<n; i++)
    for (unsigned int j=0; j<n; j++)
      for (unsigned int k=0; k<n; k++)
	{
	  // Note: Access the 1D arrays from [1] ... [n]
	  //conical_x[gp] = jacobi_x[j+1];                     
	  //conical_y[gp] = gauss_x[i+1] * (1.-jacobi_x[j+1]); 
	  //conical_w[gp] = gauss_w[i+1] * jacobi_w[j+1];      
	  
	  conical_x[gp] = jacB1D.qp(k)(0);                                                  //t[k];
	  conical_y[gp] = jacA1D.qp(j)(0)  * (1.-jacB1D.qp(k)(0));                         //s[j]*(1.-t[k]);
	  conical_z[gp] = gauss1D.qp(i)(0) * (1.-jacA1D.qp(j)(0)) * (1.-jacB1D.qp(k)(0)); //r[i]*(1.-s[j])*(1.-t[k]);
	  conical_w[gp] = gauss1D.w(i)     * jacA1D.w(j)          * jacB1D.w(k);          //A[i]*B[j]*C[k];
	  gp++;
	}
  
  // Print out the conical product points and weights in a form we can use
  // (These are now in 0-based storage)
  for (unsigned int i=0; i<n_total; ++i)
    {
      std::cout << "_points[" << std::setw(2) << i << "](0)=" << conical_x[i] << "L;\n";
      std::cout << "_points[" << std::setw(2) << i << "](1)=" << conical_y[i] << "L;\n";
      std::cout << "_points[" << std::setw(2) << i << "](2)=" << conical_z[i] << "L;\n";
    }

  std::cout << "\n";
  for (unsigned int i=0; i<n_total; ++i)
    std::cout << "_weights[" << std::setw(2) << i << "]=" << conical_w[i] << "L;\n";

  // As a check of the method, compute the sum of the weights
  mpfr_class sumweights=0.0;
  for (unsigned int i=0; i<n_total; ++i)
    sumweights += conical_w[i];

  std::cout << "Sum of weights = " << sumweights << std::endl;
  
  return 0;
}
