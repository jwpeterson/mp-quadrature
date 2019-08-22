// Functions for computing exact integrals of monomials
// on the reference triangle.
#include "exact.h"
#include "matrix.h"
#include "vect.h"
#include "newton.h"
#include "solver_data.h"

// C++ includes
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm> // std::sort
#include <numeric> // std::accumulate
#include <stdlib.h> // random
#include <time.h> // time()
#include <assert.h>

// Function used to (possibly) compute the residual and Jacobian at
// the input u.  If either "r" or "jac" is nullptr, their computation
// is skipped.
void residual_and_jacobian(std::vector<mpfr_class> * r,
                           Matrix<mpfr_class> * jac,
                           const std::vector<mpfr_class> & u);

// Use nonlinear conjugate gradient method to search for a minimum
// https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
bool nlcg(std::vector<mpfr_class> & u);

// Implements the gradient descent method for the scalar function
// f = (1/2) \vec{r} \cdot \vec{r},
// where \vec{r} is the residual vector.
bool gradient_descent(std::vector<mpfr_class> & u);

// Generates a fixed-size, dense block of random samples using the
// Latin hypercube sampling method and checks whether any of them
// converge.
void generate_lhc_and_test();

// Our goal here is to determine whether there might be more than one
// solution to the SEVENTH-order Ro3-invariant quadrature rule
// originally reported in Gatermann, 1988.  We are going to try
// several different sets of starting data and see if any of them
// converge to a different solution than the one reported in the
// paper.
//
// For more information, see gatermann.C
int main()
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Varies the sequence of pseudorandom numbers returned by random().
  // srandom(1566449865);

  // Use the current time since epoch as a seed.
  time_t seed = time(nullptr);
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  unsigned int testset_counter = 0;
  while (true)
    {
      ++testset_counter;
      std::cout << "Running test set " << testset_counter << std::endl;
      generate_lhc_and_test();
      break;
    }

  return 0;
}



void generate_lhc_and_test()
{
  // Set parameters to be used by solvers
  SolverData solver_data;
  solver_data.verbose = true;
  solver_data.residual_and_jacobian = residual_and_jacobian;

  // A "random" initial guess. Note that all parameters must be
  // positive, the sum of the weights cannot exceed the reference element
  // area. Also, each weight parameter corresponds to 3 quadrature
  // points, so actually it cannot exceed (1/3) * (1/2) = 1/6.
  // For the points, they must be chosen such that 1 - x - y > 0
  // and their relative order matters, i.e. reversing x and y will
  // give you a different solution.
  std::vector<mpfr_class> w_guess(4);
  std::vector<mpfr_class> x_guess(4);
  std::vector<mpfr_class> y_guess(4);

  // Latin hypercube sampling.
  // The total number of random numbers we need is n_dim = 11.  If the
  // total number of tests is n_tests, then we will need n_tests *
  // n_dim random numbers. We generate all the points at the same time
  // to ensure that they don't overlap.  1 / n_tests = "bin size", we
  // guarantee only one sample from each bin size over the entire set
  // of tests. Note also: n_tests = n_bins, so we only need one variable
  // for this. Note that we can't make the number of tests *too* big, since
  // we are building a dense matrix of random numbers here. So it may be
  // necessary to wrap this in an outer loop.
  unsigned int n_dim = 11;
  unsigned int n_tests = 1;
  double bin_size = 1. / double(n_tests);

  // Status message
  std::cout << "Generating random initial guesses..." << std::endl;

  // random_numbers will be n_dim * n_tests in size. Each *column*
  // will therefore represent an n_dim-dimensional set of test
  // parameters. We will random_shuffle() the rows so that we get one
  // sample from each "bin" of a given random variable.
  std::vector<std::vector<double>> random_numbers(n_dim);
  for (unsigned int i=0; i<random_numbers.size(); ++i)
    {
      random_numbers[i].resize(n_tests);

      // Get random number in the current bin.
      for (unsigned int j=0; j<random_numbers[i].size(); ++j)
        {
          // Lower bound of current bin
          double lower = j * bin_size;

          random_numbers[i][j] = bin_size * double(random()) / RAND_MAX + lower;
        }
    }

  // random_shuffle all rows after the first one. This way we don't have
  // a sample from the *same* bin for each random variable.
  for (unsigned int i=0; i<random_numbers.size(); ++i)
    std::random_shuffle(random_numbers[i].begin(), random_numbers[i].end());

  // Print with fewer digits
  // auto flags = std::cout.flags();
  // std::cout.precision(8);
  // for (unsigned int i=0; i<random_numbers.size(); ++i)
  //   {
  //     for (unsigned int j=0; j<random_numbers[i].size(); ++j)
  //       std::cout << random_numbers[i][j] << " ";
  //     std::cout << std::endl;
  //   }
  // std::cout.flags(flags);

  // Status message
  std::cout << "Testing initial guesses..." << std::endl;

  // Counter to keep track of the number of initial guesses tested.
  unsigned int guess = 0;
  unsigned int n_converged = 0;
  unsigned int n_converged_last_message = 1;

  for (unsigned int t=0; t<n_tests; ++t)
    {
      ++guess;

      if (guess == n_tests/2)
        std::cout << "Half way done..." << std::endl;

      if (n_converged > n_converged_last_message && (n_converged % 10 == 0))
        {
          std::cout << "Found " << n_converged << " duplicate converged solutions." << std::endl;
          n_converged_last_message = n_converged;
        }

      std::vector<double> current_random;
      current_random.reserve(n_dim);
      for (unsigned int j=0; j<n_dim; ++j)
        current_random.push_back(random_numbers[j][t]);

      // Debugging: for a small n_tests, make sure that each random
      // number comes from a different bin.
      // std::cout << "Random numbers used for current test." << std::endl;
      // print(current_random);

      // Every time we need a random number, post-increment "ctr".
      unsigned int ctr = 0;

      for (unsigned int i=0; i<w_guess.size(); ++i)
        {
          w_guess[i] = mpfr_class(1) / 6;

          for (unsigned int j=0; j<i; ++j)
            w_guess[i] -= w_guess[j];

          // The last weight doesn't get randomized, it's whatever is left
          // after choosing the first 3.
          if (i < w_guess.size() - 1)
            w_guess[i] *= current_random[ctr++];
        }

      // Set random initial guesses for x and y:
      for (auto & val : x_guess)
        val = current_random[ctr++];
      for (unsigned int i=0; i<y_guess.size(); ++i)
        y_guess[i] = (1. - x_guess[i]) * current_random[ctr++];

      // We are done, so make sure we used all the random numbers.
      if (ctr != n_dim)
        {
          std::cout << "Incorrect number of random number used!" << std::endl;
          abort();
        }

      // For the true solution, the sum of weight parameters is:
      // 2.65e-02 + 4.38e-02 + 2.87e-02 + 6.74e-02 = 0.1664,
      // so we should verify the same is true for our random initial
      // guess.
      // print(w_guess);
      // std::cout << "sum = "
      //           << std::accumulate(w_guess.begin(), w_guess.end(), mpfr_class(0))
      //           << std::endl;

      // Verify initial guess is OK.
      for (unsigned int i=0; i<w_guess.size(); ++i)
        {
          // std::cout << "(x" << i << ", y" << i << ")="
          //           << "(" << x_guess[i] << ", " << y_guess[i] << ")" << std::endl;

          mpfr_class barycentric_point = mpfr_class(1) - x_guess[i] - y_guess[i];
          if (barycentric_point < 0)
            {
              std::cout << "Error, invalid initial guess!" << std::endl;
              std::abort();
            }
        }

      // Set initial guess vector
      std::vector<mpfr_class> & u = solver_data.u;
      u.clear();
      u.reserve(12);
      for (unsigned int i=0; i<w_guess.size(); ++i)
        {
          u.push_back(w_guess[i]);
          u.push_back(x_guess[i]);
          u.push_back(y_guess[i]);
        }

      // Debugging: set initial guess to known solution!
      u =
        {
          2.65e-02, // w1
          6.23e-02, // x1
          6.75e-02, // y1

          4.38e-02, // w2
          5.52e-02, // x2
          3.21e-01, // y2

          2.87e-02, // w3
          3.43e-02, // x3
          6.60e-01, // y3

          6.74e-02, // w4
          5.15e-01, // x4
          2.77e-01, // y4
        };

      // Print initial guess
      // std::cout << "Initial guess=" << std::endl;
      // print(u);

      // Status flag returned by different methods.
      bool converged = false;

      // Do some number Newton minimization iterations.  If you follow
      // 25 iterations of this up by Newton-Raphson iterations, on
      // seed 1566449865 it finds 810 converged (duplicate) solutions,
      // 388 at the half-way point, on a run of 100000.  It only found
      // ~60 solutions on such a test when doing Newton-Raphson
      // iterations alone. This is a much slower approach since
      // approximating the Hessian is expensive, but it might be worth
      // it if it is much more likely to converge. Using fewer Newton
      // minimization iterations (10) after 10k iterations with the
      // same seed we found about 32 converged (duplicate)
      // solutions. This extrapolates to about 320 on 10k runs. So it
      // seems that doing more iterations does have a benefit.  We get
      // two (2) converged solutions on a test with 1k runs.  Using 25
      // minimization iteration on the same seed with n_tests=1000
      // gave four converged solutions. So doing more 2.5x as many
      // newton_min iterations gave about 2x more successful runs.
      // Unfortunately after "fixing" the Hessian I got a "Matrix is
      // singular" error at some point during the 100k run that had
      // worked fine with the wrong Hessian! By comparison, the same
      // run had 93 (half-way) and 187 (out of 100k) converged
      // solutions when doing 25 iterations of gradient_descent()
      // first. gradient_descent() does seem to run much faster than
      // newton_min().
      // Note: after fixing the Hessian, the rate of finding solutions
      // was actually worse for the same seed, only 180 were found by
      // the half-way point. This is about 2x as many as gradient descent,
      // but I'm pretty sure it's more than 2x as slow so it may not be
      // worthwhile overall.

      // Debugging: print u before minimization
      // std::cout << "u=" << std::endl;
      // print(u);

      converged = newton_min(solver_data);

      // Debugging: print u after minimization
      // std::cout << "u=" << std::endl;
      // print(u);

      // Use the gradient descent method. This is much cheaper and simpler
      // to understand than newton_min, and it seems to do about as good a job?
      // converged = gradient_descent(u);

      // Nonlinear conjugate gradient method
      // For reference, I once let this algorithm run for about 4.8M iterations
      // and it didn't reduce the residual all that much (this is a squared residual
      // so the "real" value is O(1.e-4)
      // Iteration 4780800, residual = 4.2396738648405063868090553697656e-8
      // converged = nlcg(u);

      // Follow up initial guess minimization with Newton. For a
      // random initial guess, it might be worthwhile seeing if we can
      // improve it slightly with an optimization algorithm before
      // switching to Newton's method?
      converged = newton(solver_data);

      if (converged)
        {
          // Debugging
          std::cout << "Solution converged!" << std::endl;

          // Print converged solution
          std::cout << "u=" << std::endl;
          print(u);

          // We want to check whether this solution is a permutation
          // of the one that we already have. So we sort the weights
          // and compare them to the weights of the solution we
          // already have.
          std::vector<mpfr_class> weights = {u[0], u[3], u[6], u[9]};
          std::sort(weights.begin(), weights.end());

          const std::vector<mpfr_class> known_weights =
            {
              mpfr_class("2.6517028157436251428754180460739e-2"),
              mpfr_class("2.8775042784981585738445496900219e-2"),
              mpfr_class("4.3881408714446055036769903139288e-2"),
              mpfr_class("6.7493187009802774462697086166421e-2")
            };

          for (unsigned int i=0; i<weights.size(); ++i)
            weights[i] -= known_weights[i];

          mpfr_class norm_w = norm(weights);
          // std::cout << "norm_w=" << norm_w << std::endl;

          if (norm_w > mpfr_class(1.e-6))
            {
              std::cout << "Possible new solution u=" << std::endl;
              print(u);
            }
          else
            {
              // Keep track of the number of converged (but duplicate) solutions.
              n_converged++;
            }
        }
    } // end loop over n_tests
}



bool nlcg(std::vector<mpfr_class> & u)
{
  // nlcg parameters.
  const mpfr_class tol(1.e-36);
  const mpfr_class divtol(1.e16);
  // const mpfr_class alphamin(1.e-10); // not currently used
  unsigned iter = 0;
  const unsigned int maxits = 10;
  bool converged = false;

  // Problem size
  unsigned int n = u.size();

  // Storage needed for algorithm.
  std::vector<mpfr_class> r(n), du(n), /*u_old(n),*/ du_old(n), grad_f(n), s(n), s_old(n), trial_u(n);
  Matrix<mpfr_class> jac(n,n);

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual and Jacobian at the current guess.
      residual_and_jacobian(&r, &jac, u);

      // Compute the size of the scalar function "f"
      // which we are trying to minimize, _not_ dot(r,r)^0.5.
      mpfr_class residual = 0.5 * dot(r,r);

      // Debugging:
      // if (iter % 100 == 0)
      //   std::cout << "Iteration " << iter << ", residual = " << residual << std::endl;

      // If the residual is small enough, return true.
      if (residual < tol)
        {
          converged = true;
          break;
        }

      // If the residual is too large, just give up.
      if (residual > divtol)
        break;

      // Compute the gradient of the objective function.
      grad_f = jac.matvec_transpose(r);

      // du = -1 * grad_f; // does not work for some reason?
      for (unsigned int i=0; i<n; ++i)
        du[i] = -grad_f[i];

      if (iter == 1)
        {
          // Initially we take a gradient descent step.
          s = du;
        }
      if (iter > 1)
        {
          // Debugging
          // std::cout << "du_old=" << std::endl;
          // print(du_old);
          // std::cout << "du=" << std::endl;
          // print(du);

          // Compute beta using "FR" formula
          mpfr_class beta = dot(du, du) / dot(du_old, du_old);

          // Compute beta using "PR" formula
          // mpfr_class beta = dot(du, du - du_old) / dot(du_old, du_old);

          // Compute beta using "HS" formula
          // mpfr_class beta = -dot(du, du - du_old) / dot(s_old, du - du_old);

          // Debugging:
          // std::cout << "beta = " << beta << std::endl;

          // If beta is negative then we have to set beta to zero and use
          // the s = -grad(f). The "FR" method will not give a negative beta...
          // if (beta < 0)
          //   beta = 0;

          // It does not make sense to do more than n iterations of nlcg,
          // you are supposed to take a gradient descent step (beta=0) before continuing.
          if (iter % n == 0)
            {
              // std::cout << "Taking gradient descent step to 'reset' nlcg." << std::endl;
              beta = 0;
            }

          s = du + beta * s_old;
        }

      // Compute next iterate.
      // TODO: we should actually take a step u += alpha * s, where
      // alpha is chosen by line search.
      mpfr_class alpha = mpfr_class(1);

      for (unsigned int back=0; back<10; ++back)
        {
          // std::cout << "Trying step with alpha=" << alpha << std::endl;
          trial_u = u + alpha * s;
          residual_and_jacobian(&r, nullptr, u);
          mpfr_class trial_residual = 0.5 * dot(r,r);
          if (trial_residual < residual)
            break;
          else
            alpha /= 2;
        }

      // std::cout << "Accepting step with alpha=" << alpha << std::endl;

      // In gradient descent, the following "alpha" is used:
      // if (iter > 1)
      //   {
      //     std::vector<mpfr_class> z = u - u_old;
      //     std::vector<mpfr_class> y = du_old - du; // grad_f - grad_f_old;
      //     alpha = abs(dot(z, y)) / dot(y, y);
      //   }
      u += alpha * s;

      du_old = du;
      s_old = s;
      // u_old = u;
    } // end while(true)

  return converged;
}

bool gradient_descent(std::vector<mpfr_class> & u)
{
  // Gradient descent parameters.
  const mpfr_class tol(1.e-36);
  const mpfr_class divtol(1.e16);
  // const mpfr_class alphamin(1.e-10);
  // const bool do_backtracking = true; // not currently used.
  unsigned iter = 0;
  const unsigned int maxits = 25;
  bool converged = false;

  // Problem size
  unsigned int n = u.size();

  // Storage needed for algorithm.
  std::vector<mpfr_class> r(n), du(n), grad_f(n), grad_f_old(n), u_old(n);
  Matrix<mpfr_class> jac(n,n);

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual and Jacobian at the current guess.
      residual_and_jacobian(&r, &jac, u);

      // Compute the size of the scalar function "f"
      // which we are trying to minimize, _not_ dot(r,r)^0.5.
      mpfr_class residual = 0.5 * dot(r,r);

      // Debugging:
      // std::cout << "Iteration " << iter << ", residual = " << residual << std::endl;

      // If the residual is small enough, return true.
      if (residual < tol)
        {
          converged = true;
          break;
        }

      // If the residual is too large, just give up.
      if (residual > divtol)
        break;

      // Compute grad(f) = jac^T * r at the current guess.
      grad_f = jac.matvec_transpose(r);

      // Debugging
      // std::cout << "grad_f=" << std::endl;
      // print(grad_f);

      // Compute scalar coefficient gamma that tells how far down the gradient direction to
      // travel.
      mpfr_class gamma_n(1);
      if (iter > 1)
        {
          std::vector<mpfr_class> z = u - u_old;
          std::vector<mpfr_class> y = grad_f - grad_f_old;

          gamma_n = abs(dot(z, y)) / dot(y, y);

          // Debugging:
          // std::cout << "gamma_n = " << gamma_n << std::endl;
        }

      // Store old information for next iteration. Note: do this before updating u.
      grad_f_old = grad_f;
      u_old = u;

      // Compute the new iterate.
      subtract_scaled(u, gamma_n, grad_f);
    }

  return converged;
}



void
residual_and_jacobian(std::vector<mpfr_class> * r,
                      Matrix<mpfr_class> * jac,
                      const std::vector<mpfr_class> & u)
{
  // If there's nothing to do, then there's nothing to do.
  if (!r && !jac)
    return;

  // The list of monomial exponents that we are going to test.
  std::vector<std::pair<unsigned int, unsigned int>> polys =
    {
      {0,0}, {2,0}, {3,0}, {2,1}, {4,0}, {5,0},
      {4,1}, {6,0}, {5,1}, {7,0}, {6,1}, {4,2}
    };

  // The number of equations (and the number of unknowns).
  unsigned int n = polys.size();

  // Zero any previous values and allocate space for residual and Jacobian, if required.
  if (r)
    {
      r->clear();
      r->resize(n);
    }
  if (jac)
    {
      jac->clear();
      jac->resize(n, n);
    }

  // For each basis funciton i, we loop over all (w,x,y) triples in u and
  // compute the residual and Jacobian contributions, as requested.
  for (unsigned int i=0; i<n; i++)
    {
      // The powers of x and y in the monomial we are currently integrating.
      unsigned int xpower = polys[i].first;
      unsigned int ypower = polys[i].second;

      // Accumlate quadrature rule contribution for r(i).
      for (unsigned int q=0; q<u.size(); q+=3)
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
              (*jac)(i,q+1) += w * (
                xpower * pow(x, xpm1) * pow(y, ypower) +
                (mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) + pow(z, xpower) * ypower * pow(x, ypm1)) +
                pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1));

              // Derivative wrt y
              (*jac)(i,q+2) += w * (
                pow(x, xpower) * ypower * pow(y, ypm1) +
                mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) +
                (xpower * pow(y, xpm1) * pow(z, ypower) + pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1)));
            }
        }

      // Subtract off the true integral value, I(p_i)
      if (r)
        (*r)[i] -= exact_tri(xpower, ypower);
    }
}
