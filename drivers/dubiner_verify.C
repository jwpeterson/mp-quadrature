#include "common_definitions.h"
#include "dubiner.h"
#include "conical.h"

// This test driver verifies the orthogonality of the Dubiner
// polynomials, both symbolic and numeric versions.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  Dubiner dubiner;
  Conical conical;

  // Read max Dubiner polynomial degree from command line
  unsigned int max_dubiner_degree = 4;
  if (argc > 1)
    max_dubiner_degree = atoi(argv[1]);
  std::cout << "\nVerifying d=" << max_dubiner_degree << " Dubiner polynomials." << std::endl;

  // Compute integral(phi(i)*phi(j)) for each of the Dubiner polynomials.
  for (unsigned dubiner_degree=/*0*/max_dubiner_degree; dubiner_degree <= max_dubiner_degree; ++dubiner_degree)
    {
      // The number of polynomials in the Dubiner basis of degree 'dubiner_degree'
      const unsigned Np = (dubiner_degree+1)*(dubiner_degree+2)/2;

      // The points and weights from a conical product quadrature rule
      // of high enough order to integrate phi(i)**2, which has degree 2*dubiner_degree
      conical.rule2D(2*dubiner_degree);
      const std::vector<Point<mpfr_class> > & conical_rule_points = conical.get_points();
      const std::vector<mpfr_class> & conical_rule_weights = conical.get_weights();

      // Nq is the number of points in the quadrature rule
      const unsigned Nq = conical_rule_points.size();

      // Place to store the integrated mass and Laplace matrices.
      Matrix<mpfr_class>
        mass_sums(Np, Np),
        laplace_sums(Np, Np);

      for (unsigned q=0; q<Nq; ++q)
        {
          // std::cout << "\nq=" << q << std::endl;

          std::vector<mpfr_class> current_vals;
          std::vector<Point<mpfr_class> > current_derivs;

          dubiner.p(dubiner_degree,
                    /*xi=*/  conical_rule_points[q](0),
                    /*eta=*/ conical_rule_points[q](1),
                    current_vals,
                    current_derivs);

          if (current_vals.size() != Np || current_derivs.size() != Np)
            {
              std::cerr << "current_vals.size() = " << current_vals.size() << " should be Np = " << Np << std::endl;
              std::cerr << "current_derivs.size() = " << current_derivs.size() << " should be Np = " << Np << std::endl;
              std::abort();
            }

          // std::cout << "computing mass matrix entries." << std::endl;

          // Make matrices of the phi(i)*phi(j) and dphi(i)*dphi(j) values
          Matrix<mpfr_class>
            mass(Np, Np),
            laplace(Np, Np);

          // // Extra debugging:
          // for (unsigned i=0; i<Np; ++i)
          //   {
          //     std::cout << "current_derivs[" << i << "]=" << current_derivs[i](0) << ", " << current_derivs[i](1) << std::endl;
          //     // std::cout << "current_vals[" << i << "]=" << current_vals[i] << std::endl;
          //   }

          for (unsigned i=0; i<Np; ++i)
            for (unsigned j=0; j<Np; ++j)
              {
                // std::cout << i << "," << j << std::endl;
                mass(i,j) = current_vals[i]*current_vals[j];
                laplace(i,j) = current_derivs[i](0)*current_derivs[j](0) + current_derivs[i](1)*current_derivs[j](1);
              }

          // std::cout << "accumulating quadrature sums." << std::endl;

          // Accumulate the phi(i)*phi(j) and dphi(i)*dphi(j) integrals using conical product rule quadrature.
          for (unsigned i=0; i<Np; ++i)
            for (unsigned j=0; j<Np; ++j)
              {
                mass_sums(i,j) += conical_rule_weights[q] * mass(i,j);
                laplace_sums(i,j) += conical_rule_weights[q] * laplace(i,j);
              }
        }

      // Print the mass matrix
      std::cout << "\n int phi(i)*phi(j)" << std::endl;
      mass_sums.print();

      // Print an error message if an off-diagonal entry is too large
      for (unsigned i=0; i<Np; ++i)
        for (unsigned j=0; j<Np; ++j)
          {
            if ((i != j) && (abs(mass_sums(i,j)) > 1.e-30))
              {
                std::cerr << "Matrix entry " << i << "," << j << " should be zero, but is " << mass_sums(i,j) << std::endl;
                std::abort();
              }
          }

      // Print the Laplace matrix
      std::cout << "\n int dphi(i)*dphi(j)" << std::endl;
      laplace_sums.print();
    }

  return 0;
}
