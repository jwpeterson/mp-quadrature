#include "test_ro3.h"
#include "solver_data.h"
#include "lhc.h"
#include "gradient_descent.h"
#include "nlcg.h"
#include "newton.h"
#include "vect.h"

// C++ includes
#include <stdlib.h> // random()

void test_ro3(unsigned int d,
              unsigned int n_tests,
              SolverData & solver_data)
{
  // The number of unknowns
  unsigned int N = (d*d + 3*d + 6)/6;
  // The number of centroid points (0 or 1)
  unsigned int n_centroid = N % 3;
  // The number of Ro3 points
  unsigned int n_ro3 = N / 3;

  // Debugging
  std::cout << "N = " << N << std::endl;
  std::cout << "n_centroid = " << n_centroid << std::endl;
  std::cout << "n_ro3 = " << n_ro3 << std::endl;

  if (n_centroid > 1)
    {
      std::cout << "Error, can't have a rule with more than one centroid point." << std::endl;
      abort();
    }

  // A "random" initial guess. Note that all parameters must be
  // positive, the sum of the weights cannot exceed the reference element
  // area. Also, each weight parameter corresponds to 3 quadrature
  // points, so actually it cannot exceed (1/3) * (1/2) = 1/6.
  // For the points, they must be chosen such that 1 - x - y > 0
  // and their relative order matters, i.e. reversing x and y will
  // give you a different solution.
  std::vector<mpfr_class> w_guess(n_centroid + n_ro3);
  std::vector<mpfr_class> x_guess(n_ro3);
  std::vector<mpfr_class> y_guess(n_ro3);

  // The total number of random numbers we need is 1 less than the
  // number of parameters we are solving for, since we choose the
  // weights so that they sum up to the reference element area.
  unsigned int n_dim = n_centroid + (3 * n_ro3) - 1;

  // Status message
  std::cout << "Generating random initial guesses..." << std::endl;

  // Latin hypercube sampling.
  // random_numbers will be n_dim * n_tests in size. Each *column*
  // will therefore represent an n_dim-dimensional set of test
  // parameters. We will random_shuffle() the rows so that we get one
  // sample from each "bin" of a given random variable.
  // std::vector<std::vector<double>> random_numbers =
  //   lhc(n_dim, n_tests);

  // No LHC sampling, just an (n_dim x n_tests) block of random numbers.
  std::vector<std::vector<double>> random_numbers(n_dim);
  for (auto & vec : random_numbers)
    {
      vec.resize(n_tests);
      for (auto & val : vec)
        val = double(random()) / RAND_MAX;
    }

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
          std::cout << "Incorrect number of random numbers used!" << std::endl;
          abort();
        }

      // // Verify initial guess is OK.
      // for (unsigned int i=0; i<w_guess.size(); ++i)
      //   {
      //     // std::cout << "(x" << i << ", y" << i << ")="
      //     //           << "(" << x_guess[i] << ", " << y_guess[i] << ")" << std::endl;
      //
      //     mpfr_class barycentric_point = mpfr_class(1) - x_guess[i] - y_guess[i];
      //     if (barycentric_point < 0)
      //       {
      //         std::cout << "Error, invalid initial guess!" << std::endl;
      //         std::abort();
      //       }
      //   }

      // Set initial guess vector
      std::vector<mpfr_class> & u = solver_data.u;
      u.clear();
      u.reserve(n_dim + 1);

      // If there's a centroid point, it's weight appears first in the
      // vector of unknowns.
      if (n_centroid)
        u.push_back(w_guess[0]);

      for (unsigned int i=0; i<x_guess.size(); ++i)
        {
          u.push_back(w_guess[i + n_centroid]);
          u.push_back(x_guess[i]);
          u.push_back(y_guess[i]);
        }

      // For d==3, there is a four-point rule, but it has a negative
      // centroid weight, -27/96. It is actually better to use a 2x2
      // conical product rule which has all positive weights in this
      // case.
//      if (d==3)
//        u =
//          {
//            -2.8125000000000000000000000000000e-1
//            2.6041666666666666666666666666667e-1
//            6.0000000000000000000000000000000e-1
//            2.0000000000000000000000000000000e-1
//          };

      // For d==5, we get the following 7 point rule. I believe that
      // this should be the same fifth-order rule that we currently
      // have in libmesh which is listed as being from "Quadrature on
      // Simplices of Arbitrary Dimension" by Walkington. This rule
      // actually has an analytical solution for the points and
      // weights. Note that x==y for the two Ro3 points, this
      // corresponds to points located on the three triangle medians:
      // (a,a), (a,1-2*a), and (1-2*a,a). It is also interesting to
      // note that our algorithm of using minimization/rootfinding
      // "cycles" converges for almost every random initial guess in
      // this case (although sometimes requiring more than one cycle),
      // i.e. this is a really easy problem.
//      if (d==5)
//        u =
//          {
//            1.1250000000000000000000000000000e-1,
//            6.6197076394253090368824693916576e-2,
//            4.7014206410511508977044120951345e-1,
//            4.7014206410511508977044120951345e-1,
//            6.2969590272413576297841972750091e-2,
//            1.0128650732345633880098736191512e-1,
//            1.0128650732345633880098736191512e-1
//          };

      // Debugging: set initial guess to known solution for 12 point, degree=7 case.
//      if (d==7)
//        {
//          u =
//            {
//              2.65e-02, // w1
//              6.23e-02, // x1
//              6.75e-02, // y1
//
//              4.38e-02, // w2
//              5.52e-02, // x2
//              3.21e-01, // y2
//
//              2.87e-02, // w3
//              3.43e-02, // x3
//              6.60e-01, // y3
//
//              6.74e-02, // w4
//              5.15e-01, // x4
//              2.77e-01, // y4
//            };
//
//          // Check that if we change the order of the orbits, this
//          // is still a solution...
//          std::cout << "u before permutation=" << std::endl;
//          print(u);
//
//          // Permute orbits 1 and 2
//          Matrix<mpfr_class> P12(u.size(), u.size());
//          P12.get_values() =
//            {
//              0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
//              1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
//            };
//
//          // Permute orbits 1 and 3
//          Matrix<mpfr_class> P13(u.size(), u.size());
//          P13.get_values() =
//            {
//              0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
//              0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
//              1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
//              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
//            };
//
//          // Apply a sequence of permutations
//          u = P12 * u;
//          u = P13 * u;
//
//          std::cout << "u after permutation=" << std::endl;
//          print(u);
//        }

      // This rule is reported in Cools & Haegemans TW96 for degree=8
      // but it seems quite dubious to me as one of the weights is
      // basically zero and one of the points is waaaay outside the
      // reference element... it does converge using our solver
      // however, and the solution is:
      //
      // u=
      // 1.6058343856681218798400231921758e-10
      // 3.4579201116826902882140270303267e-1
      // 3.6231682215692616666791285362864
      // 2.6530624434780379346899858317533e-2
      // 6.5101993458939166328104791860010e-2
      // 8.7016510156356306077747992432900e-1
      // 2.9285717640165892159260292342121e-2
      // 6.5177530364879570753723544017162e-1
      // 3.1347788752373300717357808540986e-1
      // 4.3909556791220782401863560411239e-2
      // 3.1325121067172530695595743477398e-1
      // 6.3062143431895614010295743491499e-1
      // 6.6940767639916174191830767611771e-2
      // 5.1334692063945414949358896494129e-1
      // 2.8104124731511039057273791035622e-1
// if (d==8)
//      u =
//        {
//          0.16058e-09, // w1
//          0.3458, // x1
//          0.3623e+01, // y1
//
//          0.2653e-01, // w2
//          0.651e-01, // x2
//          0.8701, // y2
//
//          0.2928e-01, // w3
//          0.6518, // x3
//          0.3135, // y3
//
//          0.4391e-01, // w4
//          0.31325, // x4
//          0.6306, // y4
//
//          0.6694e-01, // w5
//          0.513347, // x5
//          0.2810, // y5
//        };


      // This rule is reported in Cools & Haegemans TW96 for degree=10,
      // but it also has points outside the reference element. Our solver
      // converges to this if we do Newton iterations only, but not if we
      // do some minimization steps first..
      //
      // u=
      // 4.7910534861520060666258795267317e-2
      // 1.5319130036758557630224218647182e-7
      // 5.8469201683584513030840399896305e-2
      // -5.4887778772527519316973900501509e-1
      // 1.3260526227928785221366811484040e-2
      // 5.0849285064031410704623553266007e-2
      // 9.0799059794957813439289694951544e-1
      // 1.5646439344539042136182802356194e-2
      // 5.1586732419949674486519820945603e-1
      // 4.6312452842927062902198120980987e-1
      // 2.1704258224807323310774404349324e-2
      // 2.4311033191739048229513038265738e-1
      // 7.2180595182371959467025855850715e-1
      // 2.1797613600129922367395409006103e-2
      // 7.5397765920922660134372366560338e-1
      // 2.0647569839132397632996966382170e-1
      // 3.8587913508193459468029657331332e-2
      // 4.2209207910846960293600275784774e-1
      // 1.2689533413411127326858250178468e-1
      // 3.9699584282594413021921681475048e-2
      // 1.9823878846663354067849545395741e-1
      // 6.2124412566393319744542137613089e-1
// if (d==10)
//      u =
//        {
//          0.4791e-01, // w0
//
//          0.15319e-06, // w1
//          0.5846e-01, // x1
//          -0.5489, // y1
//
//          0.1326e-01, // w2
//          0.50849e-01, // x2
//          0.908, // y2
//
//          0.15646e-01, // w3
//          0.515867e+00, // x3
//          0.46312, // y3
//
//          0.2170e-01, // w4
//          0.24311e+00, // x4
//          0.7218, // y4
//
//          0.218e-01, // w5
//          0.75398e+00, // x5
//          0.2065, // y5
//
//          0.3859e-01, // w6
//          0.422e+00, // x6
//          0.1269, // y6
//
//          0.397e-01, // w7
//          0.1982e+00, // x7
//          0.6212, // y7
//        };

      // Print initial guess
      // std::cout << "Initial guess=" << std::endl;
      // print(u);

      // Status flag returned by different methods.
      bool converged = false;

      // Try a few iterations of one of the minimization routines.
      // converged = newton_min(solver_data);
      // converged = gradient_descent(solver_data);
      // converged = nlcg(solver_data);

      // Follow up initial guess minimization with Newton. For a
      // random initial guess, it might be worthwhile seeing if we can
      // improve it slightly with an optimization algorithm before
      // switching to Newton's method?
      // converged = newton(solver_data);

      // Cycle back and forth between N minimization steps and N Newton iterations.
      // while (true)
      unsigned int n_cycles = 1;
      mpfr_class norm_r = 0.;
      for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
        {
          // Minimization step
          // gradient_descent(solver_data);
          // newton_min(solver_data);
          // nlcg(solver_data);

          // Root-finding step
          converged = newton(solver_data);

          // Report one residual value per cycle
          if (n_cycles > 1)
            {
              std::vector<mpfr_class> r;
              solver_data.residual_and_jacobian(&r, nullptr, u);
              norm_r = norm(r);
              std::cout << "cycle " << cycle << ", residual norm = " << norm_r << std::endl;
            }
          // print(r);

          if (converged)
            break;
        }

      if (converged)
        {
          // Print converged solution
          std::cout << "Solution converged!" << std::endl;
          std::cout << "u=" << std::endl;
          print(u);

          // Keep track of the number of converged solutions.
          n_converged++;
        }

      // Debugging: we didn't converge, but maybe it was still promising?
      if (!converged)
        {
          std::cout << "Solution is *not* converged!" << std::endl;
          std::cout << "u=" << std::endl;
          print(u);
        }

    } // end loop over n_tests
}
