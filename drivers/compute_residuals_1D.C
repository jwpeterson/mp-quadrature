#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory> // std::auto_ptr
#include <getopt.h> // getopt_long()
#include "legendre.h"
#include "gauss.h"

// Function which computes and prints ||r_N||_2 values for the 1D rule
// defined by the points x and weights w.  The 'filebase' string is used
// in naming the output line plot files.  Returns the final H1-norm of the
// error (but note that it prints the H1-norm squared while it is running!)
mpfr_class compute_rn_norms(const std::vector<mpfr_class> & x,
                            const std::vector<mpfr_class> & w,
                            const std::string & filebase);


// Class representing a function, its exact integral, L2, and H1 norms.
// We assume that the integrand has at least one parameter named x0.
struct Integrand
{
  // Constructor
  Integrand(const mpfr_class & x0_in) :
    x0(x0_in)
  {}

  // Destructor
  virtual ~Integrand() {}

  // The parameter
  mpfr_class x0;

  // Returns the value of f at the point x
  virtual mpfr_class f(const mpfr_class & x) const = 0;

  // Returns the exact integral f
  virtual mpfr_class exact() const = 0;

  // Returns the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const = 0;

  // Returns the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const = 0;

  // Returns the H1-norm.
  mpfr_class H1_norm() const
  {
    return sqrt(this->L2_norm_squared() + this->H1_semi_norm_squared());
  }
};



// Class representing the function |x-x_0|^{alpha}. The function is in
// H1 if alpha>1/2.
struct AlphaSingularity : Integrand
{
  AlphaSingularity(const mpfr_class & x0_in,
                   const mpfr_class & alpha_in) :
    Integrand(x0_in),
    alpha(alpha_in)
  {
    if (alpha <= mpfr_class(1.)/2.)
      {
        std::cerr << "alpha must be greater than 1/2." << std::endl;
        std::abort();
      }
  }

  mpfr_class alpha;

  // f = abs(x-x_0)^{\alpha}, \alpha > \frac{1}{2}
  virtual mpfr_class f(const mpfr_class & x) const
  {
    return pow(abs(x - x0), alpha);
  }

  // I(f) = 1/(\alpha+1) [ (1+x_0)^{\alpha+1} + (1-x_0)^{\alpha+1} ]
  virtual mpfr_class exact() const
  {
    return mpfr_class(1.)/(alpha+1.) * (pow(1+x0, alpha+1) + pow(1-x0, alpha+1));
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return mpfr_class(1.)/(2.*alpha+1.) * (pow(1+x0, 2.*alpha+1) + pow(1-x0, 2.*alpha+1));
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    return alpha*alpha/(2.*alpha - 1.) * (pow(1+x0, 2.*alpha-1) + pow(1-x0, 2.*alpha-1));
  }
};



// Class representing the piecewise integrand
// { 1 + (x-x0)/(1+x0)
// { 1 - (x-x0)/(1-x0)
struct Piecewise : Integrand
{
  Piecewise(const mpfr_class & x0_in) :
    Integrand(x0_in)
  {}

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const
  {
    // Handle edge cases without dividing by zero.
    if (x0 == -1.)
      return 0.5*(1.-x);
    if (x0 == 1.)
      return 0.5*(1.+x);

    // Handle in-between cases
    if (x <= x0)
      return mpfr_class(1.) + (x-x0)/(1.+x0);
    else
      return mpfr_class(1.) - (x-x0)/(1.-x0);
  }

  // I(f) = 1
  virtual mpfr_class exact() const
  {
    return 1.;
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return mpfr_class(2.)/3.;
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    // Handle edge cases without dividing by zero.
    if (x0 == -1. || x0 == 1.)
      return 0.5;

    return mpfr_class(1.)/(1. + x0) + mpfr_class(1.)/(1. - x0);
  }
};




// Class representing the exponential integrand
// exp(-alpha*|x-x0|)
struct Exponential : Integrand
{
  Exponential(const mpfr_class & x0_in,
              const mpfr_class & alpha_in) :
    Integrand(x0_in),
    alpha(alpha_in)
  {}

  mpfr_class alpha;

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const
  {
    return exp(-alpha*abs(x-x0));
  }

  // True integral
  virtual mpfr_class exact() const
  {
    return 2./alpha * (1. - exp(-alpha)*cosh(alpha*x0));
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return 1./alpha * (1. - exp(-2.*alpha)*cosh(2.*alpha*x0));
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    return alpha * (1. - exp(-2.*alpha)*cosh(2.*alpha*x0));
  }
};



// Class representing the sinusoidal integrand
// sin(alpha*(x-x0))
struct Sinusoidal : Integrand
{
  Sinusoidal(const mpfr_class & x0_in,
             const mpfr_class & alpha_in) :
    Integrand(x0_in),
    alpha(alpha_in)
  {}

  mpfr_class alpha;

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const
  {
    return sin(alpha*(x-x0));
  }

  // True integral
  virtual mpfr_class exact() const
  {
    return -2./alpha * sin(alpha) * sin(alpha * x0);
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return 1. - sin(2*alpha) * cos(2*alpha*x0) / (2.*alpha);
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    return alpha*alpha * (1. + sin(2*alpha) * cos(2*alpha*x0) / (2.*alpha));
  }
};





// Approximate the value of an integral on [-1,1] using the quadrature
// rule defined by the points x, and weights w.  Returns the absolute value
// of the quadrature error, |E(f)|.
mpfr_class test_integral(const std::vector<mpfr_class> & x,
                         const std::vector<mpfr_class> & w,
                         const Integrand & integrand);


// A usage function for this utility
void usage();


// Sets up the rule to be used with e.g. compute_rn_norms() or test_integral()
void set_up_rule(unsigned rule_number,
                 std::vector<mpfr_class> & x,
                 std::vector<mpfr_class> & w,
                 std::string & filebase);



int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Default rule number
  unsigned rule_number = 0;

  // options descriptor - this associates "long" and "short" command line options.
  static struct option longopts[] =
    {
      {"rule-number", required_argument, NULL, 'r'},
      {"help",        no_argument,       NULL, 'h'},
      { NULL,         0,                 NULL,  0 }
    };

  // Parse command line options using getopt_long()
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "hr:", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'r':
          {
            rule_number = atoi(optarg);
            break;
          }

        case 'h':
          {
            usage();
            return 0;
          }

        default:
          // We could error here, print a usage command, or just ignore unrecognized arguments...
          usage();
        }
    } // end while

  // Set up the quadrature rule to use
  std::vector<mpfr_class> x, w;
  std::string filebase;
  set_up_rule(rule_number, x, w, filebase);

  // Compute the ||r_N||_1 values
  mpfr_class norm_rN = compute_rn_norms(x, w, filebase);
  std::cout << "Final ||r_N||_1 value = " << norm_rN << std::endl;

  // Set up the integrand, we have a couple to choose from...
  std::auto_ptr<Integrand> integrand;

  // To build up the name of the output file to write the abs_err and err_bound values to
  // TODO: if the plots/ directory does not exist, we must create it first, otherwise nothing
  // will be written.
  std::ostringstream oss;

  unsigned integrand_type = 3;
  switch (integrand_type)
    {
    case 0:
      {
        mpfr_class alpha = mpfr_class(2.)/3.;
        integrand.reset(new AlphaSingularity(/*x0=*/0., alpha));
        oss << "plots/" << filebase << "_alpha_" << std::setprecision(5) << alpha << "_err_bounds.csv";
        break;
      }

    case 1:
      {
        integrand.reset(new Piecewise(/*x0=*/0.));
        oss << "plots/" << filebase << "_piecewise_err_bounds.csv";
        break;
      }

    case 2:
      {
        integrand.reset(new Exponential(/*x0=*/0., /*alpha=*/6.));
        oss << "plots/" << filebase << "_exponential_err_bounds.csv";
        break;
      }

    case 3:
      {
        // If alpha=n*pi, then the integrand does not depend on x0.
        // const_pi() is defined in gmpfrxx.h.  The choice alpha=3*pi
        // is particularly interesting: it causes all 3 of the
        // first-order quadrature rules to give *exactly* the same
        // result.  The Gauss(1) and Closed Trap. rules always have
        // exactly the same (abs) error when alpha=n*pi, although the
        // closed-trap one flips back and forth in sign.  This is
        // because the true integrand is zero when alpha=n*pi.
        mpfr_class alpha = 2. * const_pi();

        // The integrand I(f) is maximized (for x0 = +/-1) when alpha=atan(2)
        // mpfr_class alpha = atan(mpfr_class(2.));

        // A large number which is not an integer multiple of pi
        // mpfr_class alpha = 100.;

        integrand.reset(new Sinusoidal(/*x0=*/0., alpha));
        oss << "plots/" << filebase << "_sinusoidal_alpha_" << std::setprecision(5) << alpha << "_err_bounds.csv";
        break;
      }

    default:
      std::abort();
    }

  // Open file. We don't need crazy precision on the err bounds.
  std::ofstream out(oss.str().c_str());
  out.precision(6);
  out.setf(std::ios_base::scientific);

  // Compute the abs_err and err_bound values at several plot points
  unsigned n_plot_points = 200;
  mpfr_class dx = mpfr_class(2.) / mpfr_class(n_plot_points-1);
  for (unsigned i=0; i<n_plot_points; ++i)
    {
      // Move the singularity to the next point
      integrand->x0 = -1. + i*dx;

      mpfr_class abs_err = test_integral(x, w, *integrand);

      // Compute the product ||f||_1 * ||r_N||_1.  Our a priori estimate
      // says |E(f)| must be less than this (plus some other terms we hope
      // are small).
      mpfr_class norm_f = integrand->H1_norm();

      if (abs_err > norm_f*norm_rN)
        {
          std::cerr << "Oops, the actual error was larger than the a prior bound!" << std::endl;
          std::abort();
        }

      // Print x0, abs_err, ||f||_1, ||r_N||_1 (note: last column is a constant)
      out << integrand->x0 << ", " << abs_err << ", " << norm_f << ", " << norm_rN << ",\n";
    }

  return 0;
}



void set_up_rule(unsigned rule_number,
                 std::vector<mpfr_class> & x,
                 std::vector<mpfr_class> & w,
                 std::string & filebase)
{
  switch (rule_number)
    {
    case 0:
      {
        // The (closed) trapezoidal rule is exact for linears
        x.resize(2); w.resize(2);
        x[0] = -1;   w[0] = 1.;
        x[1] =  1;   w[1] = 1.;
        filebase = "trapezoidal_1st_closed";
        break;
      }

    case 1:
      {
        // The 1-point Gauss rule is exact for linears
        gauss_rule(/*n=*/1, x, w);
        filebase = "gauss_1st";
        break;
      }

    case 2:
      {
        // The (open) trapezoidal rule is also exact for linears
        x.resize(2); w.resize(2);
        x[0] = mpfr_class(-1.)/3.;   w[0] = 1.;
        x[1] = mpfr_class( 1.)/3.;   w[1] = 1.;
        filebase = "trapezoidal_1st_open";
        break;
      }

      // --------------------------------------------------------------------------------
      // 3rd-order rules
      // --------------------------------------------------------------------------------

    case 3:
      {
        // Milne(3) rule
        // A 3rd-order accurate rule with 3-points, basically the open NC
        // equivlent of Simpson's rule.  It has a negative weight.
        x.resize(3); w.resize(3);
        x[0] = -0.5;  w[0] = mpfr_class(4.)/3.;
        x[1] =  0.0;  w[1] = mpfr_class(-2.)/3.;
        x[2] =  0.5;  w[2] = mpfr_class(4.)/3.;
        filebase = "milne_3rd";
        break;
      }

    case 4:
      {
        // Points/weights for the Simpson(3) rule
        x.resize(3); w.resize(3);
        x[0] = -1;   w[0] = mpfr_class(1.)/3.;
        x[1] =  0;   w[1] = mpfr_class(4.)/3.;
        x[2] =  1;   w[2] = mpfr_class(1.)/3.;
        filebase = "simpson_3rd";
        break;
      }

    case 5:
      {
        // Gauss(2) rule
        gauss_rule(/*n=*/2, x, w);
        filebase = "gauss_3rd";
        break;
      }

    case 6:
      {
        // Open N-C(4) rule
        x.resize(4); w.resize(4);
        x[0] = mpfr_class(-3.)/5;  w[0] = mpfr_class(11.)/12.;
        x[1] = mpfr_class(-1.)/5.; w[1] = mpfr_class(1. )/12.;
        x[2] = mpfr_class( 1.)/5.; w[2] = mpfr_class(1. )/12.;
        x[3] = mpfr_class( 3.)/5.; w[3] = mpfr_class(11.)/12.;
        filebase = "noname_3rd";
        break;
      }

    case 7:
      {
        // Simpson(4) rule
        // AKA Simpson's "3/8" rule.  It has 4 points, but is still
        // 3rd-order, so presumably no one uses it because of this.  Note
        // that the name comes from the factor of h/8 out in front, since
        // our h==2, it's more like the 3/4 rule.
        x.resize(4); w.resize(4);
        x[0] = -1;                 w[0] = mpfr_class(1.)/4.;
        x[1] = mpfr_class(-1.)/3.; w[1] = mpfr_class(3.)/4.;
        x[2] = mpfr_class( 1.)/3.; w[2] = mpfr_class(3.)/4.;
        x[3] =  1;                 w[3] = mpfr_class(1.)/4.;
        filebase = "simpson_3rd_three_eighths";
        break;
      }

      // --------------------------------------------------------------------------------
      // 5th-order rules
      // --------------------------------------------------------------------------------

    case 8:
      {
        // Open-NC(5) rule
        x.resize(5); w.resize(5);
        x[0] = mpfr_class(-2.)/3.; w[0] = mpfr_class( 11.)/10.;
        x[1] = mpfr_class(-1.)/3.; w[1] = mpfr_class(-14.)/10.;
        x[2] = 0.0;                w[2] = mpfr_class( 26.)/10.;
        x[3] = mpfr_class( 1.)/3.; w[3] = mpfr_class(-14.)/10.;
        x[4] = mpfr_class( 2.)/3.; w[4] = mpfr_class( 11.)/10.;
        filebase = "opennc_5th";
        break;
      }

    case 9:
      {
        // Gauss(3) rule
        gauss_rule(/*n=*/3, x, w);
        filebase = "gauss_5th";
        break;
      }

    case 10:
      {
        // 5-point "Boole's rule", exact for quintics
        x.resize(5); w.resize(5);
        x[0] = -1.0;   w[0] = mpfr_class(7. )/45.;
        x[1] = -0.5;   w[1] = mpfr_class(32.)/45.;
        x[2] =  0.0;   w[2] = mpfr_class(12.)/45.;
        x[3] =  0.5;   w[3] = mpfr_class(32.)/45.;
        x[4] =  1.0;   w[4] = mpfr_class(7. )/45.;
        filebase = "boole_5th";
        break;
      }

    default:
      {
        std::cerr << "Unrecognized rule_number = " << rule_number << std::endl;
        std::abort();
      }
    }
}



void usage()
{
  std::cout << "\n";
  std::cout << "This program computes discrete residual norm values and tests the quadrature rule with a representative integrand.\n";
  std::cout << "\n";
  std::cout << "Valid command line options are:\n";
  std::cout << "--help, -h\n";
  std::cout << "  Display this help message\n";
  std::cout << "--rule-number, -r N\n";
  std::cout << "  Test rule N, where:\n";
  std::cout << "    Third-order rules:\n";
  std::cout << "      0: Closed Trapezoidal (2)\n";
  std::cout << "      1: Gauss (1)\n";
  std::cout << "      2: Open Trapezoidal (2)\n";
  std::cout << "    Third-order rules:\n";
  std::cout << "      3: Milne (3)\n";
  std::cout << "      4: Simpson (3)\n";
  std::cout << "      5: Gauss (2)\n";
  std::cout << "      6: Open N-C (4)\n";
  std::cout << "      7: Simpson (4)\n";
  std::cout << "    Fifth-order rules:\n";
  std::cout << "      8: Open N-C (5)\n";
  std::cout << "      9: Gauss (3)\n";
  std::cout << "     10: Boole (5)\n";
  std::cout << "\n";
}



mpfr_class test_integral(const std::vector<mpfr_class> & x,
                         const std::vector<mpfr_class> & w,
                         const Integrand & integrand)
{
  mpfr_class sum = 0.;
  for (unsigned q=0; q<x.size(); ++q)
    sum += w[q] * integrand.f(x[q]);

  // std::cout << "sum=" << sum << std::endl;

  // Compute the exact integrand value
  mpfr_class exact = integrand.exact();
  // std::cout << "exact=" << exact << std::endl;

  // Compute norms
  // std::cout << "L2-norm=" << sqrt(integrand.L2_norm_squared()) << std::endl;
  // std::cout << "H1-seminorm=" << sqrt(integrand.H1_semi_norm_squared()) << std::endl;
  // std::cout << "L2-norm**2=" << integrand.L2_norm_squared() << std::endl;
  // std::cout << "H1-seminorm**2=" << integrand.H1_semi_norm_squared() << std::endl;
  // std::cout << "H1-norm=" << integrand.H1_norm() << std::endl;

  // Compute abs(E(f)) by subtracting
  return abs(sum-exact);
}


mpfr_class compute_rn_norms(const std::vector<mpfr_class> & x,
                            const std::vector<mpfr_class> & w,
                            const std::string & filebase)
{

  // Verify integration of Dubiner polynomials.  Note: they should all
  // be zero except the first one, which should be 1.0
  Legendre legendre;

  // The degree of Legendre polynomials to compute the residual with
  const unsigned max_legendre_degree = 30;

  // Eventually we'll return the square root of this
  mpfr_class residual_h1_norm_squared = 0.;

  for (unsigned legendre_degree=0; legendre_degree <= max_legendre_degree; ++legendre_degree)
    {
      std::vector<mpfr_class> E, current_vals;

      for (unsigned i=0; i<x.size(); ++i)
        {
          // Evaluate all the Legendre polynomials at the current qp
          legendre.p(legendre_degree, x[i], current_vals);

          // After the first loop iteration, this resize() should do nothing
          E.resize(current_vals.size());

          for (unsigned j=0; j<current_vals.size(); ++j)
            E[j] += w[i] * current_vals[j];
        }

      // E[] now contains all of the E(phi_i) except the first entry,
      // which needs to have the "volume" of the region, 2,
      // subtracted.  This is because the true values of all the
      // integrals is zero except for the first one.
      E[0] -= 2.;

      // Debugging: print E[i] values
      //  std::cout << "E=" << std::endl;
      //  for (unsigned i=0; i<E.size(); ++i)
      //    std::cout << E[i] << std::endl;

      // Build the H1-projection matrix
      Matrix<mpfr_class> A;
      legendre.build_H1_projection_matrix(legendre_degree, A);

      // Debugging: print A, verify that it is symmetric, etc.  Note that
      // A.print(); isn't very good at lining up the columns, so we print
      // by hand to verify...
      if (false)
        {
          // save the initial format flags
          std::ios_base::fmtflags flags = std::cout.flags();

          std::cout.precision(8);

          for (unsigned i=0; i<A.n_rows(); ++i)
            {
              for (unsigned j=0; j<A.n_rows(); ++j)
                {
                  // Convert to double for more control over printing...
                  // the mpfr_class stuff ignores formatting for '0' for
                  // some reason.
                  std::cout << A(i,j).get_d() << " ";
                }
              std::cout << std::endl;
            }

          // restore the original format flags
          std::cout.flags(flags);
        }

      // Make sure there are the same number of coeffs as residual entries
      if (E.size() != A.n_rows())
        {
          std::cerr << "E vector size does not match matrix size!" << std::endl;
          std::abort();
        }

      // The LU solve will destroy A, therefore make a copy
      // of A first since we need it to compute the H1-norm of r.
      Matrix<mpfr_class> A_copy(A);
      // std::cout<< "A_copy.n_rows() = " << A_copy.n_rows() << std::endl;

      // Solve Ar = E for r.
      std::vector<mpfr_class> r;
      A.lu_solve(r, E);

      // Compute the H1-norm (squared) of the residual.  To do this, we need
      // the orthogonality coeffs, or in general something like A^2
      residual_h1_norm_squared = 0.;
      for (unsigned i=0; i<r.size(); ++i)
        for (unsigned j=0; j<r.size(); ++j)
          residual_h1_norm_squared += r[i] * A_copy(i,j) * r[j];

      std::cout << fix_string(residual_h1_norm_squared, /*append_L=*/false) << std::endl;

      // For non-zero residuals, we can also plot the r_N function by
      // computing sum r[i] * phi[i] at a bunch of
      // linearly-spaced points...
      if (residual_h1_norm_squared > 1.e-30)
        {
          std::ostringstream oss;
          oss << "plots/" << filebase << "_line_plot_" << std::setw(2) << std::setfill('0') << legendre_degree;
          std::ofstream out(oss.str().c_str());

          unsigned n_plot_points = 200;
          mpfr_class dx = mpfr_class(2.) / mpfr_class(n_plot_points-1);
          for (unsigned i=0; i<n_plot_points; ++i)
            {
              mpfr_class current_x = -1. + i*dx;

              // Evaluate the basis functions at current_x
              legendre.p(legendre_degree, current_x, current_vals);

              // Sum
              mpfr_class sum = 0;
              for (unsigned j=0; j<current_vals.size(); ++j)
                sum += r[j] * current_vals[j];

              // Print comma-separated x and y values for line plot
              out << fix_string(current_x, /*append_L=*/false) << ", " << fix_string(sum, /*append_L=*/false) << std::endl;
            }
        }
    }

  return sqrt(residual_h1_norm_squared);
}
