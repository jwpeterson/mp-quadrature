// System headers
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdlib.h> // random()
#include <getopt.h> // getopt_long()

// The GNU multi-precision library
#include "gmpfrxx.h"

// Header files for this project
#include "common_definitions.h"
#include "matrix.h"

// Print the usage information for this driver program.
void usage();

// This driver program implements LU decomposition and back
// substitution using multi-precision variables.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // Define what type of scalars we will use.  The initial goal
  // is to keep the code working with built-in types as well, this
  // is a good way to test that it does still compile and run...
  typedef mpfr_class scalar_type;
  // typedef Real scalar_type;

  // options descriptor - this associates several "long" and "short"
  // command line options.  The last element of the longopts array has
  // to be filled with zeros.
  static struct option longopts[] =
    {
      {"matrix-size",     required_argument, NULL, 'n'},
      {"binary-digits",   required_argument, NULL, 'b'},
      {"help",            no_argument,       NULL, 'h'},
      { NULL,             0,                 NULL,  0 }
    };

  // Problem size
  unsigned int n = 3;

  // Number of binary digits to use
  // 53 binary digits is about what you normally get with a double.
  unsigned n_binary_digits = 256;

  // Parse command line options using getopt_long()
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "xs:ub:r:e:io:ph", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'n':
          n = atoi(optarg);
          break;

        case 'b':
          n_binary_digits = atoi(optarg);
          break;

        case 'h':
          usage();
          return 0;

        default:
          usage();
        }
    } // end while

  // Set default mpfr_class precision
  mpfr_set_default_prec(n_binary_digits);

  // Row-major matrix
  Matrix<scalar_type> A(n,n);

  // manufactured solution, solution, and rhs vectors
  std::vector<scalar_type> manufactured_x(n), x(n), b(n);

  // Examples of ill-conditioned matrices
  // .) Vandermonde
  // .) Octave reports cond ~ 10^5 for this little matrix:
  //    [-149, -50, -154; 537, 180, 546; -27, -9, -25]
  // .) A = 1/2*[1 1; 1+1.e-10 1-1.e-10 ], b=[1;1], x=[1;1]
  //    Octave reports cond ~ 10^10 for this matrix
  // .) Hilbert matrix (http://en.wikipedia.org/wiki/Hilbert_matrix)
  //    H(i,j) = 1/(i+j-1), i,j=1,2,...n
  //    The 5x5 Hilbert matrix has cond ~ 10^5
  //    In Octave, hilb(N) generates the N*N Hilbert matrix.
  //    octave-3.4.0:17> cond(hilb(6))
  //    ans =  1.49510586407410e+07
  //    octave-3.4.0:18> cond(hilb(7))
  //    ans =  4.75367356437253e+08
  //    octave-3.4.0:19> cond(hilb(8))
  //    ans =  1.52575757800826e+10

  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      {
        // Fill matrix with random [-1,1] values
        // A(i,j) = (2.*(static_cast<Real>(random()) / static_cast<Real>(RAND_MAX)) - 1.);

        // Fill matrix with Hilbert matrix entries
        A(i,j) = static_cast<scalar_type>(1.) / static_cast<scalar_type>( (i+1) + (j+1) - 1 );
      }


  // Print A (is this printing garbage or actual random digits in the case of mpfr_class?)
  if (n <= 5)
    {
      std::cout << "A=" << std::endl;
      std::cout << A;
    }

  // Fill up manufactured_x with values, compute A*manufactured_x = b, then attempt LU solve
  // with A, b, verify that you can get back the original manufactured_x.
  for (unsigned int i=0; i<n; ++i)
    {
      // Random values in [-1,1]
      manufactured_x[i] = (2.*(static_cast<Real>(random()) / static_cast<Real>(RAND_MAX)) - 1.);
    }

  // Print manufactured_x
  std::cout << "\nmanufactured_x=" << std::endl;
  for (unsigned int i=0; i<n; ++i)
    {
      std::cout << manufactured_x[i] << std::endl;
    }

  // Compute A*manufactured_x, store result in b
  for (unsigned int i=0; i<n; ++i)
    {
      b[i] = 0.0;
      for (unsigned int j=0; j<n; ++j)
        b[i] += A(i,j)*manufactured_x[j];
    }

  // Print b
  std::cout << "\nb=" << std::endl;
  for (unsigned int i=0; i<n; ++i)
    std::cout << b[i] << std::endl;

  // Solve, using LU, the sytem Ax=b for x
  A.lu_solve(x,b);

  if (n <= 5)
    {
      std::cout << "After LU solve: " << std::endl;
      std::cout << A;
    }

  std::cout << "\nx=" << std::endl;
  for (unsigned int i=0; i<n; ++i)
    std::cout << x[i] << std::endl;

  // Compare x and manufactured_x
  scalar_type l2_norm = 0.;
  for (unsigned int i=0; i<n; ++i)
    {
      l2_norm += (x[i]-manufactured_x[i])*(x[i]-manufactured_x[i]);
    }

  std::cout << "\n||x-manufactured_x||, Discrete L2-norm =" << sqrt(l2_norm) << std::endl;

  // Results for Hilbert example, discrete error in ||x-manufactured_x||
  // n   MPFR,512   MPFR,256  double
  // 6   2.72e-148  8.13e-70  6.92e-12
  // 10  3.42e-142  5.86e-66  1.29e-4
  // 12  1.98e-139  1.78e-62  3.85e-1
  // 14  1.14e-136  1.99e-59  1.25e+01

  // This can be confired with a simple Matlab/Octave script as well...
  // A=hilb(10);
  // manufactured_x=rand(10,1);
  // b=A*manufactured_x;
  // x = A\b;
  // norm(x-manufactured_x,2)
  // ans =  9.2225e-04

  return 0;
}



void usage()
{
  std::cout << "\n";
  std::cout << "This program tests the LU solve capability of our Matrix class.\n";
  std::cout << "\n";
  std::cout << "Valid command line options are:\n";
  std::cout << "--matrix-size, -n #     = Solve a matrix of size n*n.\n";
  std::cout << "--binary-digits, -b #   = Default number of binary digits to use for mpfr_class variables.\n";
  std::cout << "--help, -h              = Print this message.\n";
  std::cout << "\n";
}
