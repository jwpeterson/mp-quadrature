// nlopt includes. The preprocessor define is set on the compile line
// with a -D flag.
#ifndef HAVE_NLOPT

int main()
{
  std::cout << "This driver requires nlopt. Install nlopt and then rerun configure, "
            << "using the --with-nlopt-include and --with-nlopt-lib flags to specify "
            << "its location." << std::endl;
  return 0;
}

#else

// mp-quadrature includes
#include "solver_data.h"
#include "newton.h"
#include "gradient_descent.h"
#include "nlcg.h"
#include "ro3.h"

// nlopt includes
#include <nlopt.h>

// C/C++ includes
#include <iostream>
#include <math.h>
#include <stdlib.h> // random()
#include <getopt.h> // getopt_long()

// Objective function f : R^n -> R and, if grad != nullptr, grad(f).
double myfunc(unsigned n, const double * x, double * grad, void * my_func_data);

// Constraints are of the form fc(x) <= 0, where fc is required to have
// the same parameters as the objective function. In our case, there is
// one inequality constraint associated with each ro3 orbit of the form:
// 1 - x - y >= 0,
// which we multiply by -1 to write as:
// x + y - 1 <= 0
// So in this case, the incoming "data" pointer just needs to tell us
// the index of the x dof for the orbit being constrained.
double myconstraint(unsigned n, const double *x, double *grad, void *data);

// Convert a numerical nlopt_result successful return code to a human readable string.
std::string nlopt_result_to_string(nlopt_result res);

// Call this function when the program is run with the wrong arguments.
void usage();

int main(int argc, char ** argv)
{
  // You can't trust all these digits from doubles, but you can with
  // mpfr_class objects.
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // Process command line arguments
  static struct option longopts[] =
    {
      {"degree",      required_argument, NULL, 'd'},
      {"n-centroid",  required_argument, NULL, 'c'},
      {"n-vertex",    required_argument, NULL, 'v'},
      {"n-edge",      required_argument, NULL, 'e'},
      {"n-general",   required_argument, NULL, 'g'},
      {"help",        no_argument,       NULL, 'h'},
      { NULL,         0,                 NULL,  0 }
    };

  // To be set from the command line
  unsigned int d = 0;

  // Number of centroid, vertex, edge, and general orbits.
  unsigned int nc = 0;
  unsigned int nv = 0;
  unsigned int ne = 0;
  unsigned int ng = 0;

  // A colon following an option means it has an argument.
  // If there's an unrecognized argument, getopt_long()
  // prints a message, so we don't handle it in the cases below.
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "hd:c:v:e:g:", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'd': { d = atoi(optarg); break; }
        case 'c': { nc = atoi(optarg); break; }
        case 'v': { nv = atoi(optarg); break; }
        case 'e': { ne = atoi(optarg); break; }
        case 'g': { ng = atoi(optarg); break; }
        case 'h': { usage(); return 0; }
        }
    } // end while

  if (d == 0 || nc > 1 || nv > 1)
    {
      std::cout << "Error, invalid command line arguments provided!" << std::endl;
      usage();
      return 1;
    }

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  // It's faster to comment this out while debugging.
  mpfr_set_default_prec(256);

  // Use the current time since epoch as a seed.
  time_t seed = time(nullptr);
  // long seed = 1567251384;
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  // Build ro3 rule object
  Ro3 r(d, nc, nv, ne, ng);

  // d==2, dim==2
  // Note: This case is somewhat tricky despite also being the
  // simplest non-trivial case!  The issue is that the Jacobian is
  // singular at the root, so Newton iterations converge very slowly,
  // and you don't obtain the full requested accuracy in the final
  // solution.  I am a bit surprised that lu_solve() does not throw a
  // singular matrix exception, I think the problem is that the matrix
  // is only approximately singular, and we only throw if it's
  // *exactly* singular.
  // -d2 -c0 -v0 -e1 -g0 # 3 QP <-- Solution confirmed

  // d==3, dim==4
  // -d3 -c1 -v0 -e0 -g1 # 4 QP <-- No solution (-ve wt soln only)
  // -d3 -c0 -v0 -e2 -g0 # 6 QP <-- No solution
  // -d3 -c1 -v1 -e1 -g0 # 7 QP <-- New (?) solution

  // d==4, dim==5
  // o The best known deg=4 PI rule is a D3-invariant rule with 6 QP.
  // o We also found a degree=4 rule with 6 QPs having ne=ng=1.
  // o For the case with nc=1, ne=2 there does not seem to be a solution, but
  //   the same (global) minimum objective function value of ~2.57337-08
  // is found by the algorithm regardless of starting point.
  // -d4 -c0 -v0 -e1 -g1 # 6 QP <-- New (?) soln found!
  // -d4 -c1 -v1 -e0 -g1 # 7 QP <-- No solution
  // -d4 -c1 -v0 -e2 -g0 # 7 QP <-- No solution
  // -d4 -c0 -v1 -e2 -g0 # 9 QP <-- No solution

  // d==5, dim==7
  // -d5 -c1 -v0 -e0 -g2 # 7 QP <-- Also a D3 rule in libmesh

  // d=6, dim==10
  // -d6 -c1 -v0 -e0 -g3 # 10 QP <-- No solution
  // -d6 -c0 -v1 -e0 -g3 # 12 QP <-- No solution
  // -d6 -c0 -v0 -e2 -g2 # 12 QP <-- New (?) solution
  // -d6 -c1 -v1 -e1 -g2 # 13 QP <-- No solution
  // -d6 -c1 -v0 -e3 -g1 # 13 QP <-- No solution
  // -d6 -c0 -v0 -e5 -g0 # 15 QP <-- No solution
  // -d6 -c1 -v1 -e4 -g0 # 16 QP <-- No solution

  // d==7, dim=12
  // -d7 -c0 -v0 -e0 -g4 # 12 QP <-- Gatermann solution
  // -d7 -c1 -v0 -e1 -g3 # 13 QP <-- New (?) solution
  // -d7 -c0 -v1 -e1 -g3 # 15 QP <-- No solution
  // -d7 -c0 -v0 -e3 -g2 # 15 QP <-- New (?) solution
  // -d7 -c1 -v1 -e2 -g2 # 16 QP <-- New (?) solution
  // -d7 -c1 -v0 -e4 -g1 # 16 QP <-- No solution
  // -d7 -c0 -v1 -e4 -g1 # 18 QP <-- No solution
  // -d7 -c1 -v1 -e5 -g0 # 18 QP <-- No solution

  // d==8, dim=15, best PI degree 8 rule in libmesh has 16 QPs
  // -d8 -c0 -v0 -e0 -g5 # 15 QP <-- No solution
  // -d8 -c1 -v0 -e1 -g4 # 16 QP <-- New (?) solution
  // -d8 -c0 -v1 -e1 -g4 # 18 QP <-- New (?) solution
  // -d8 -c0 -v0 -e3 -g3 # 18 QP <-- No solution
  // -d8 -c1 -v1 -e2 -g3 # 19 QP <-- New (?) solution
  // -d8 -c1 -v0 -e4 -g2 # 19 QP <-- No solution
  // -d8 -c0 -v1 -e4 -g2 # 21 QP <-- No solution
  // -d8 -c0 -v0 -e6 -g1 # 21 QP <-- No solution
  // -d8 -c1 -v1 -e5 -g1 # 22 QP <-- No solution
  // -d8 -c1 -v0 -e7 -g0 # 22 QP <-- No solution
  // -d8 -c0 -v1 -e7 -g0 # 24 QP <-- No solution

  // d==9, dim=19, best PI degree 9 rule in libmesh has 19 QPs
  // -d9 -c1 -v0 -e0 -g6 # 19 QP <-- No solution
  // -d9 -c0 -v1 -e0 -g6 # 21 QP <-- New (?) solution
  // -d9 -c0 -v0 -e2 -g5 # 21 QP <-- No solution
  // -d9 -c1 -v1 -e1 -g5 # 22 QP <-- New (?) solution
  // -d9 -c1 -v0 -e3 -g4 # 22 QP <-- No solution
  // -d9 -c0 -v1 -e3 -g4 # 24 QP <-- New (?) solution
  // -d9 -c0 -v0 -e5 -g3 # 24 QP <-- No solution
  // -d9 -c1 -v1 -e4 -g3 # 25 QP <-- No solution
  // -d9 -c1 -v0 -e6 -g2 # 25 QP <-- No solution
  // -d9 -c0 -v1 -e6 -g2 # 27 QP <-- No solution
  // -d9 -c0 -v0 -e8 -g1 # 27 QP <-- No solution
  // -d9 -c1 -v1 -e7 -g1 # 28 QP <-- No solution
  // -d9 -c1 -v0 -e9 -g0 # 28 QP <-- No solution
  // -d9 -c0 -v1 -e9 -g0 # 30 QP <-- No solution

  // d==10, dim=22, best PI degree 10 rule in libmesh has 25 QPs
  // -d10 -c1 -v0 -e0 -g7 # 22 QP <-- No solution
  // -d10 -c0 -v1 -e0 -g7 # 24 QP <-- No solution
  // -d10 -c0 -v0 -e2 -g6 # 24 QP <-- TWO New (?) solutions
  // -d10 -c1 -v1 -e1 -g6 # 25 QP <-- New (?) solution
  // -d10 -c1 -v0 -e3 -g5 # 25 QP <-- New (?) solution
  // -d10 -c0 -v1 -e3 -g5 # 27 QP <-- New (?) solution
  // -d10 -c0 -v0 -e5 -g4 # 27 QP <-- No solution
  // -d10 -c1 -v1 -e4 -g4 # 28 QP <-- 4.9e-17
  // -d10 -c1 -v0 -e6 -g3 # 28 QP
  // -d10 -c1 -v0 -e7 -g2 # 28 QP
  // -d10 -c0 -v1 -e7 -g2 # 30 QP
  // -d10 -c0 -v0 -e8 -g2 # 30 QP
  // -d10 -c1 -v0 -e9 -g1 # 31 QP
  // -d10 -c0 -v1 -e9 -g1 # 33 QP
  // -d10 -c1 -v1 -e10 -g0 # 34 QP

  // d==11, dim=26, best PI degree 11 rule in libmesh has 30 QPs
  // -d11 -c0 -v0 -e1 -g8 # 27 QP <-- 1.36e-21, instance-3
  // -d11 -c1 -v1 -e0 -g8 # 28 QP <-- 7.47e-23, instance-4
  // -d11 -c1 -v0 -e2 -g7 # 28 QP <-- 1.66e-20, instance-5
  // -d11 -c0 -v1 -e2 -g7 # 30 QP <-- 3.05e-21, instance-6
  // -d11 -c0 -v0 -e4 -g6 # 30 QP <-- 3.05e-19, instance-7
  // -d11 -c1 -v1 -e3 -g6 # 31 QP <-- 4.23e-21, instance-8
  // -d11 -c1 -v0 -e5 -g5 # 31 QP <-- 4.21e-20, instance-1
  // -d11 -c0 -v1 -e5 -g5 # 33 QP <-- 3.97e-19
  // -d11 -c0 -v0 -e7 -g4 # 33 QP <-- 9.54e-16
  // -d11 -c1 -v1 -e6 -g4 # 34 QP <-- 3.23e-16
  // -d11 -c1 -v0 -e8 -g3 # 34 QP
  // -d11 -c0 -v1 -e8 -g3 # 36 QP
  // -d11 -c0 -v0 -e10 -g2 # 36 QP
  // -d11 -c1 -v1 -e9  -g2 # 37 QP
  // -d11 -c1 -v0 -e11 -g1 # 37 QP
  // -d11 -c0 -v1 -e11 -g1 # 39 QP
  // -d11 -c0 -v0 -e13 -g0 # 39 QP
  // -d11 -c1 -v1 -e12 -g0 # 40 QP

  // d==12, dim=31, best PI degree 12 rule in libmesh has 33 QPs
  // -d12 -c1 -v0 -e0 -g10 # 31 QP <-- instance-2
  // -d12 -c0 -v1 -e0 -g10 # 33 QP
  // -d12 -c0 -v0 -e2 -g9 # 33 QP
  // -d12 -c1 -v1 -e1 -g9 # 34 QP
  // -d12 -c1 -v0 -e3 -g8 # 34 QP
  // -d12 -c0 -v1 -e3 -g8 # 36 QP
  // -d12 -c0 -v0 -e5 -g7 # 36 QP
  // ...

  // d==13, dim=35, best PI degree 13 rule in libmesh has 37 QPs
  // -d13 -c0 -v0 -e1 -g11 # 36 QP
  // -d13 -c1 -v1 -e0 -g11 # 37 QP
  // -d13 -c1 -v0 -e2 -g10 # 37 QP
  // -d13 -c0 -v1 -e2 -g10 # 39 QP
  // -d13 -c0 -v0 -e4 -g9 # 39 QP
  // -d13 -c1 -v1 -e3 -g9 # 40 QP
  // -d13 -c1 -v0 -e5 -g8 # 40 QP
  // -d13 -c0 -v1 -e5 -g8 # 42 QP
  // ...

  // d==14, dim=40, best PI degree 14 rule in libmesh has 42 QPs
  // -d14 -c1 -v0 -e0 -g13 # 40 QP
  // -d14 -c0 -v0 -e2 -g12 # 42 QP
  // -d14 -c1 -v1 -e1 -g12 # 43 QP
  // -d14 -c1 -v0 -e3 -g11 # 43 QP
  // -d14 -c0 -v1 -e3 -g11 # 45 QP
  // -d14 -c0 -v0 -e5 -g10 # 45 QP
  // -d14 -c1 -v1 -e4 -g10 # 46 QP
  // ...

  // d==18, dim=64, there is no PI degree 18 rule in libmesh, next highest has 73 pts.
  // -d18 -c1 -v0 -e0 -g21 # 64 QP
  // -d18 -c0 -v1 -e0 -g21 # 66 QP
  // -d18 -c0 -v0 -e2 -g20 # 66 QP
  // -d18 -c1 -v1 -e1 -g20 # 67 QP
  // -d18 -c1 -v0 -e3 -g19 # 67 QP
  // -d18 -c0 -v1 -e3 -g19 # 69 QP
  // -d18 -c0 -v0 -e5 -g18 # 69 QP
  // -d18 -c1 -v1 -e4 -g18 # 70 QP
  // -d18 -c1 -v0 -e6 -g17 # 70 QP
  // -d18 -c0 -v1 -e6 -g17 # 72 QP
  // -d18 -c0 -v0 -e8 -g16 # 72 QP
  // -d18 -c1 -v1 -e7 -g16 # 73 QP
  // -d18 -c1 -v0 -e9 -g15 # 73 QP
  // -d18 -c0 -v1 -e9 -g15 # 75 QP
  // ...

  // Print information.
  // std::cout << "Rule has " << r.n_qp() << " quadrature points." << std::endl;
  //
  // std::cout << "CENTROID: ["
  //           << r.begin(Ro3::CENTROID)
  //           << ","
  //           << r.end(Ro3::CENTROID)
  //           << ")" << std::endl;
  //
  // std::cout << "VERTEX: ["
  //           << r.begin(Ro3::VERTEX)
  //           << ","
  //           << r.end(Ro3::VERTEX)
  //           << ")" << std::endl;
  //
  // std::cout << "EDGE: ["
  //           << r.begin(Ro3::EDGE)
  //           << ","
  //           << r.end(Ro3::EDGE)
  //           << ")" << std::endl;
  //
  // std::cout << "GENERAL: ["
  //           << r.begin(Ro3::GENERAL)
  //           << ","
  //           << r.end(Ro3::GENERAL)
  //           << ")" << std::endl;
  //
  // return 0;

  SolverData solver_data(r);
  solver_data.verbose = true;
  solver_data.maxits = 50;
  solver_data.do_backtracking = true;
  solver_data.residual_reduction_required = true;

  // Tolerances used by the optimization routines (and the 'inner'
  // optimization routine, if using AUGLAG). Recall that these are
  // objective function values, so they are proportional to the "true"
  // residual squared, i.e. an ftol of 1.e-16 corresponds to a true
  // residual of about 1.e-8. If maxeval is negative, then this
  // criterion is disabled.
  //
  // To give some idea about the relationship between maxeval and time
  // spent optimizing, we ran a few seeds of a large (-d18) test case
  // and recorded how long they took. In each case, the 'convergence'
  // reason was NLOPT_MAXEVAL_REACHED. I am not sure whether
  // specifying maxeval on the "local" optimizer effectively squares
  // the number of evaluations that are done, or if the value is
  // somehow cumulative.  Note that there seems to be a lot of
  // variability possible even for a single value of maxeval. In fact,
  // there can be overlap in the times required for 5k and 10k
  // evals. Clearly if you make the maxeval too small you will not
  // obtain any starting points close to anything like the global
  // minimum. I think the main point here is to use a reasonably
  // large value which is not infinite, so you don't spend too much
  // time in local minima where you can't achieve the required ftol.
  //
  // maxeval   Time(s) for "nlopt_test -d18 -c1 -v0 -e0 -g21"
  // 10000  &  83.72 & 36.37 & 45.09 & 81.59 & 38.72 & 142.19
  //  5000  &  19.86 & 47.02 & 19.28 & 41.47 & 18.37
  //  2500  &  9.60  & 10.32 & 9.58  & 9.46  & 21.34
  //  1000  &  3.85  & 9.23  & 9.49  & 9.37  & 6.01
  double ftol_rel = 1.e-16;
  double xtol_rel = 1.e-8;
  int maxeval = 10000;

  // "GN" = global derivative-free optimization
  // nlopt_algorithm alg = NLOPT_GN_DIRECT_L;

  // "LD" = local gradient based optimization
  // nlopt_algorithm alg = NLOPT_LD_MMA; // method of moving asymptotes
  // nlopt_algorithm alg = NLOPT_LD_SLSQP; // sequential quadratic programming
  // nlopt_algorithm alg = NLOPT_LD_TNEWTON_PRECOND_RESTART; // Preconditioned truncated Newton
  // nlopt_algorithm alg = NLOPT_LD_LBFGS; // Low storage BFGS

  // You must use a local/subsidiary optimization algorithm with AUGLAG,
  // this is set by calling nlopt_set_local_optimizer().
  nlopt_algorithm alg = NLOPT_AUGLAG;

  // MLSL also requires a local optimizer. MLSL is a multistart
  // algorithm that performs a sequence of local optimizations, and
  // uses a clustering heuristic to avoid repeated searches of the
  // same local optima. Unfortunately, I let it run for a long time
  // (over night) and it never produced any output even on the d==7
  // problem which has a known solution, so I'm not really sure what
  // it's doing or if it's even working?
  // nlopt_algorithm alg = NLOPT_G_MLSL;

  // The problem dimension depends only on "d".
  unsigned int dim = r.dim();

  // Set minimization algorithm and dimensionality
  nlopt_opt opt = nlopt_create(alg, dim);

  if (alg == NLOPT_AUGLAG ||
      alg == NLOPT_G_MLSL)
    {
      // The objective function, bounds, and nonlinear-constraint
      // parameters of local_opt are ignored. A copy is made when
      // calling set_local_optimizer, so it is safe to let the local
      // copy be destroyed.
      nlopt_opt local_opt = nlopt_create(NLOPT_LD_SLSQP, dim);
      nlopt_set_xtol_rel(local_opt, xtol_rel);
      nlopt_set_ftol_rel(local_opt, ftol_rel);
      nlopt_set_maxeval(local_opt, maxeval);
      nlopt_set_local_optimizer(opt, local_opt);
    }

  // lower bounds: all parameters must be positive, so we just use the
  // default-constructed values of the std vector for lb. All weights
  // for 3-point orbits must be less than 1/6, which is (1/3)*(ref. element area).
  // Everything else must be less than 1.
  std::vector<double> lb;
  std::vector<double> ub;
  r.bounds(lb, ub);

  // Debugging: Print entries of ub
  // std::cout << "ub=" << std::endl;
  // for (const auto & val : ub)
  //   std::cout << val << std::endl;

  nlopt_set_lower_bounds(opt, lb.data());
  nlopt_set_upper_bounds(opt, ub.data());
  nlopt_set_min_objective(opt, myfunc, static_cast<void *>(&solver_data));

  // Add inequality constraint. We need call the "add" function once
  // per inequality constraint, and each time pass it the index of the
  // x DOF of the orbit which is to be constrained.
  std::vector<unsigned int> data;
  r.inequality_constraint_indices(data);

  // Debugging:
  // std::cout << "data=" << std::endl;
  // print(data);

  for (unsigned int c=0; c<data.size(); ++c)
    nlopt_add_inequality_constraint(opt, myconstraint, &data[c], 1e-8);

  // Set tolerance(s)
  nlopt_set_xtol_rel(opt, xtol_rel);
  nlopt_set_ftol_rel(opt, ftol_rel);

  // You can use this parameter to control how long some of the
  // gradient free optimization algorithms run, otherwise it seems
  // that they can just run forever.
  nlopt_set_maxeval(opt, maxeval);

  // Solution vector
  std::vector<double> x;

  while (true)
    {
      // Generate random initial guess.
      r.guess(x);

      // // Debugging: Print entries of x
      // std::cout << "initial guess=" << std::endl;
      // for (const auto & val : x)
      //   std::cout << val << std::endl;

      // One can also solve for this rule "by hand" since if you
      // assume the edge orbit is at 1/2 it becomes a linear system of
      // equations that has a solution. This rule is not optimal since
      // it has both points on the boundary and it has 7 QPs while the
      // best known degree 3 rule has only 4 QPs. However, the weights
      // are all positive and it was a useful test case for making
      // sure that the VERTEX residual and Jacobian contributions were
      // correct.
      // if (r.d==3 && r.nc==1 && r.ne==1 && r.nv==1 && r.ng==0)
      //   x =
      //     {
      //       2.2500000000000000000000000000000e-1,
      //       2.5000000000000000000000000000000e-2,
      //       6.6666666666666666666666666666667e-2,
      //       5.0000000000000000000000000000000e-1
      //     };

      // This rule is different from the D3-invariant rule with 6 QP
      // in libmesh which is from Lyness and Jesperson and has all
      // points inside the reference element.  The rule also has 6
      // points (ne=1, ng=1) but it contains 3 points on the element
      // boundary, so the D3 rule should probably be considered superior,
      // but this one seems .
      // if (r.d==4 && r.ne==1 && r.ng==1)
      //   x =
      //     {
      //       3.4009490250701580549995389768470e-2,
      //       8.5705066495465124974799173176085e-1,
      //       1.3265717641596508611667127689820e-1,
      //       5.6577725671689651088265387031247e-1,
      //       3.1774517297791196146365253427748e-1
      //     };

      // if (r.d==5 && r.nc==1 && r.ng==2)
      //   x =
      //     {
      //       1.1250000000000000000000000000000e-1,
      //       6.6197076394253090368824693916576e-2,
      //       4.7014206410511508977044120951345e-1,
      //       5.9715871789769820459117580973105e-2,
      //       6.2969590272413576297841972750091e-2,
      //       1.0128650732345633880098736191512e-1,
      //       7.9742698535308732239802527616975e-1
      //     };

      // This rule has 12 points, which matches the degree=6 rule
      // which is currently in libmesh. There are some points on the
      // boundary of the domain, so for that reason, the rule
      // currently in libmesh may be preferred.
      // if (r.d==6 && r.ne==2 && r.ng==2)
      //   x =
      //     {
      //       1.7521923484831786607547758325659e-2,
      //       4.2749900628191084402313725002908e-1,
      //       1.2850039868925906983627292049321e-2,
      //       9.3271642599256179913317342394557e-2,
      //       5.2253795110990941510714178014166e-2,
      //       6.8836320947862802695823181103657e-2,
      //       1.9226287054542767293536118369964e-1,
      //       8.4040908201918031564777438277521e-2,
      //       5.0269934990534993982516857981038e-1,
      //       3.3904897391721909524901496402551e-1
      //     };

      // The d=7 solution with many digits of accuracy is:
      // 2.6517028157436251428754180460739e-2
      // 6.2382265094402118173683000996350e-2
      // 6.7517867073916085442557131050869e-2
      // 4.3881408714446055036769903139288e-2
      // 5.5225456656926611737479190275645e-2
      // 3.2150249385198182266630784919920e-1
      // 2.8775042784981585738445496900219e-2
      // 3.4324302945097146469630642483938e-2
      // 6.6094919618673565761198031019780e-1
      // 6.7493187009802774462697086166421e-2
      // 5.1584233435359177925746338682643e-1
      // 2.7771616697639178256958187139372e-1

      // if (r.d==7)
      //   x =
      //     {
      //       2.65e-02, // w1
      //       6.23e-02, // x1
      //       6.75e-02, // y1
      //
      //       4.38e-02, // w2
      //       5.52e-02, // x2
      //       3.21e-01, // y2
      //
      //       2.87e-02, // w3
      //       3.43e-02, // x3
      //       6.60e-01, // y3
      //
      //       6.74e-02, // w4
      //       5.15e-01, // x4
      //       2.77e-01, // y4
      //     };

      // A degree=7 rule with 13 QPs and some points on the boundary.
      // if (r.d==7 && r.nc==1 && r.ne==1 && r.ng==3)
      //   x =
      //     {
      //       5.5952307797500451168270351180298e-2,
      //       1.3472905782326158579409369641855e-2,
      //       2.9983830132009239524945782161234e-1,
      //       2.6082365542691963824789940662658e-2,
      //       5.9482416053645952551714087077197e-2,
      //       6.9397373107702483247241789818035e-2,
      //       6.1273891079961226273752784819268e-2,
      //       2.7907750657353493179879419573566e-1,
      //       1.4598341071586341030182362058933e-1,
      //       4.7186734995853834265957787816119e-2,
      //       3.2916411424668018651106653409879e-1,
      //       6.1215138856035843807465228900117e-1
      //     };

      // A degree=7 rule with 15 QPs and some points on the boundary.
      // if (r.d==7 && r.ne==3 && r.ng==2)
      //   x =
      //     {
      //       5.3055935561341320656442591841393e-3,
      //       3.4023228823368393399065798430885e-2,
      //       1.5002198243061481721648841389365e-2,
      //       7.7426490649533686580770326141160e-1,
      //       1.5301251644897155959135909968714e-2,
      //       4.0051406706170526286311293473985e-1,
      //       5.0624560660856324047249434886360e-2,
      //       1.5395946225900545445477142492632e-1,
      //       1.0329490144591944867168415690346e-1,
      //       8.0433062561717572872988221238088e-2,
      //       3.8824142455567101056261725619719e-1,
      //       4.6001634230332020766420547596794e-1
      //     };

      // A degree=7 rule with 16 QPs and points on boundary/at vertices.
      // if (r.d==7 && r.nc==1 && r.nv==1 && r.ne==2 && r.ng==2)
      //   x =
      //     {
      //       8.8156307800262736848734269834507e-2,
      //       2.0181692212805178377771076491231e-3,
      //       1.1842518154963341192378271866334e-2,
      //       8.3929917227292351714406960839597e-1,
      //       1.1842518154963341192378271866334e-2,
      //       1.6070082772707648285593039160403e-1,
      //       5.3251995907086993738381694116003e-2,
      //       6.1698507712375957239912339414237e-2,
      //       4.6915074614381202138004383029288e-1,
      //       5.8326029294951560422839897890704e-2,
      //       6.9012787955247933305015129005951e-1,
      //       1.5493606022376033347492435497024e-1
      //     };

      // A degree=8 rule with 16 QPs
      // if (r.d==8 && r.nc==1 && r.nv==0 && r.ne==1 && r.ng==4)
      //   x =
      //     {
      //       7.2521345610590588615251673993725e-2,
      //       1.1691671982623978980576722781465e-2,
      //       2.8464284047341892779748152552693e-1,
      //       1.6036995852050015675610346503080e-2,
      //       5.3867543687557559308134383474885e-2,
      //       4.6837555921867442934220084396571e-2,
      //       4.7920586242366783286185404347777e-2,
      //       4.3827550158929177203507813868026e-1,
      //       4.8025920827024682057103085428712e-1,
      //       1.4800720141846362085949314479152e-2,
      //       7.4347327164768248134793287452843e-1,
      //       1.4096246060080109128363068092795e-2,
      //       5.2042910577582663766594320557284e-2,
      //       1.8711966318476081140504944380348e-1,
      //       1.5372274637782641123563856802741e-1
      //     };

      // A degree=8 rule with 18 QPs
      // if (r.d==8 && r.nc==0 && r.nv==1 && r.ne==1 && r.ng==4)
      //   x =
      //     {
      //       1.6365619450731091079309677113646e-3,
      //       7.1563040840844654360044305010581e-3,
      //       1.4150003160281202221467222411207e-1,
      //       2.3032399676650253713746312132884e-2,
      //       2.6135943412979344118211349907110e-2,
      //       3.9683870899759689901048321440939e-1,
      //       3.0412879070999470272247013794143e-2,
      //       5.8705218537414718746555765593064e-2,
      //       1.2867350137855184316304358427429e-1,
      //       6.3077071453596363648364086313746e-2,
      //       5.0034691425209789094544109307179e-1,
      //       1.9824402124136825067504259045131e-1,
      //       4.1351450436263004488373856213471e-2,
      //       6.3013360188613761901960007277632e-1,
      //       2.9639453351548959482253455585386e-1
      //     };

      // A degree=8 rule with 19 QPs
      // if (r.d==8 && r.nc==1 && r.nv==1 && r.ne==2 && r.ng==3)
      //   x =
      //     {
      //       5.3249913011017943647121653650309e-2,
      //       1.5093537559558725922092008419750e-3,
      //       7.7922345016508128137115446990742e-3,
      //       8.5151321170633136497105403382591e-1,
      //       1.2773850763700670089813161350187e-2,
      //       4.1035393094358068585396592117871e-1,
      //       5.2019626337460031874380610292933e-2,
      //       3.0858419507892152432388479445908e-1,
      //       1.4167589821801641930724853282077e-1,
      //       2.8196727788590522799411409283777e-2,
      //       1.2245778857082282250180202616490e-1,
      //       5.6094505881456623061225804438646e-2,
      //       4.6624902515636108614766855648618e-2,
      //       6.1512488037728468907811464063495e-1,
      //       8.1799808561019703180235467829817e-2
      //     };

      // A degree=9 rule with 21 QPs
      // if (r.d==9 && r.nc==0 && r.nv==1 && r.ne==0 && r.ng==6)
      // x =
      //   {
      //     4.6208136672890810190208010916353e-4,
      //     1.3937873649405077572538903674386e-2,
      //     4.4412127969532946609530810618352e-2,
      //     6.0455484530038181374983232537504e-2,
      //     8.7959612804017479710722128519961e-3,
      //     1.9013225698280544912811185627099e-1,
      //     7.1987546384528925078927829097248e-3,
      //     4.0344560584550655466935190288604e-2,
      //     2.1110808416324593232253596761464e-1,
      //     1.2631554612265104469150870016886e-1,
      //     2.5068384777481441427534757854270e-2,
      //     4.4513058891472071846505101597941e-1,
      //     3.5535480200215007101851885821557e-2,
      //     2.5175973969753207576301584782348e-2,
      //     4.0997052576549121963973333096097e-2,
      //     2.4960151776578923653305787716817e-1,
      //     5.2881831038345628550381937105899e-2,
      //     1.9128515995430339172415385144765e-1,
      //     3.5897434142891728391398528200056e-1
      //   };

      // A degree=9 rule with 22 QPs
      // if (r.d==9 && r.nc==1 && r.nv==1 && r.ne==1 && r.ng==5)
      //   x =
      //     {
      //       5.6681242229959782378455068998313e-2,
      //       1.7709563164108618654073199527745e-4,
      //       8.0056858564321250888415067519545e-3,
      //       5.0000000000000000000000000000000e-1,
      //       2.1641769688644688644688644688645e-2,
      //       3.6838412054736283634817598783385e-2,
      //       2.2196298916076569567510252769319e-1,
      //       1.3032090619145091983675537032338e-2,
      //       4.6947443199090612267241036420043e-2,
      //       4.6947443199090612267241036420043e-2,
      //       4.1203790999884165472995247646690e-2,
      //       1.9187191273744911850366265609742e-1,
      //       1.9187191273744911850366265609742e-1,
      //       2.1641769688644688644688644688645e-2,
      //       2.2196298916076569567510252769319e-1,
      //       3.6838412054736283634817598783385e-2,
      //       4.2070716772288226519084664197013e-2,
      //       1.0044132362596745752372601888727e-1,
      //       4.4977933818701627123813699055636e-1
      //     };

      // A degree=9 rule with 24 QPs
      // if (r.d==9 && r.nc==0 && r.nv==1 && r.ne==3 && r.ng==4)
      //   x =
      //     {
      //       9.9920681638586186356375443283050e-4,
      //       7.0690247201981272962308384289864e-3,
      //       6.3193040073722041295255627701824e-1,
      //       5.9266081839999399947176520476464e-3,
      //       8.8043673507833069274213166314801e-1,
      //       5.3475597462147498356454535340421e-3,
      //       1.0667687351399658302468746396423e-1,
      //       3.7760534137683438155770128636229e-2,
      //       8.9980550840639517183429370342513e-2,
      //       3.2711502822800339810214903476533e-1,
      //       2.9572113473164001765589944002256e-2,
      //       1.0777246753344816261464172868087e-1,
      //       1.0651976200065317286577734289541e-1,
      //       5.2240014546826852011186560781300e-2,
      //       2.1197643637448513468424883083615e-1,
      //       4.8594896708812698156314636765236e-1,
      //       2.7751605042193695743962334803376e-2,
      //       6.2512240599886847137981590902263e-1,
      //       3.3384503238421439815486250637193e-1
      //     };

      // A degree=10 rule with 24 QPs. New minimum!
      // if (r.d==10 && r.nc==0 && r.nv==0 && r.ne==2 && r.ng==6)
      //   x =
      //     {
      //       5.3681429742859643492357377689009e-3,
      //       1.8721976304020603799148510523827e-1,
      //       2.2815303655641663532038669323421e-3,
      //       3.4075334397348650193113634396732e-2,
      //       4.6508993462059758854967647559714e-2,
      //       3.4147086142038659869812874556596e-1,
      //       2.0360915320008500090421409688573e-1,
      //       1.0919981154490145628957387886100e-2,
      //       6.9817400501731927826050298390783e-1,
      //       1.4413404224312299237271155349378e-2,
      //       2.4240746515267806539218234118407e-2,
      //       4.2633514579462347877517092365945e-1,
      //       3.6779605207993416066366804528172e-2,
      //       1.4936119213925489681656723297916e-2,
      //       4.0093527438242910369389026226690e-2,
      //       9.1103646481204342896318426601411e-2,
      //       3.4963363159025821525120057092851e-2,
      //       1.1690104736983285255079342004697e-1,
      //       2.7867154732579039483251817674802e-1,
      //       2.7447789822047513734307012010436e-2,
      //       7.2000132708127838543972050558105e-1,
      //       1.9134246734700412201903713737393e-1
      //     };

      // A *second* degree=10 rule with 24 QPs and 0026 configuration.
      // Also a new minimum!
      // if (r.d==10 && r.nc==0 && r.nv==0 && r.ne==2 && r.ng==6)
      //   x =
      //     {
      //       6.8034549450389437911542695284968e-3,
      //       2.4033695082490661641610501527861e-1,
      //       6.4701145825025661268470771647561e-3,
      //       5.2540484162127156065220999286902e-1,
      //       7.6297245866351294827778042218373e-3,
      //       1.2970322053864413659668358505098e-2,
      //       1.4081107816864061630771237989169e-1,
      //       2.7968592259403195267388456279467e-2,
      //       7.5136392619669750371252757075214e-1,
      //       1.5004696166167503172923442728480e-1,
      //       3.6263197047755420520828693175881e-2,
      //       5.3735517659978976913561285260804e-1,
      //       3.7661987230500154446477586010322e-1,
      //       7.4674601101896288064569309205665e-3,
      //       4.4107239622531995584271061559211e-2,
      //       2.6087950887797941493159400743790e-2,
      //       2.6160683863563246187991128905119e-2,
      //       6.5311292067804412162841018713658e-1,
      //       5.9875144970959509735752771006843e-2,
      //       4.7903439271578536483222306470543e-2,
      //       2.8533958139225425086186103213101e-1,
      //       4.8509072998273930259505659173245e-1
      //     };

      // A degree=10 rule with 25 QPs
      // if (r.d==10 && r.nc==1 && r.nv==1 && r.ne==1 && r.ng==6)
      //   x =
      //     {
      //       4.7309440928094852885976950813341e-2,
      //       7.4443068205713019768945426025270e-4,
      //       3.1829359254805817751202442366943e-3,
      //       8.3061025308964767122329536525606e-2,
      //       3.8799021945743061001599802530126e-2,
      //       1.7544304433743764192225628495153e-1,
      //       2.0933705012838168942077572054598e-1,
      //       2.0030412127900619224360889008975e-2,
      //       3.2176407820770162393255046010956e-2,
      //       3.0726941261349314805052750203576e-1,
      //       1.5850177031042592646271274554700e-2,
      //       3.9026767455703841640294956578970e-2,
      //       9.9070711082074272239198884743969e-2,
      //       1.3553817908665389393448328866720e-2,
      //       1.9269035323026781372640492307404e-2,
      //       5.7154160680610615243446739606247e-1,
      //       2.0804506634455471119732622715725e-2,
      //       1.9770352662977282642262230494035e-1,
      //       4.9561536726352818108215112762139e-2,
      //       3.7931550768623537013118400222361e-2,
      //       4.1147963472022346886180725166791e-1,
      //       1.2777722809830669309582013618709e-1
      //     };

      // A second degree=10 rule with 25 QPs
      // if (r.d==10 && r.nc==1 && r.nv==0 && r.ne==3 && r.ng==5)
      //   x =
      //     {
      //       4.0232634070269439301667832239894e-2,
      //       3.3422311840204575942302033783202e-3,
      //       5.2253403465953518825868105768782e-2,
      //       6.1424211433885610965819016355015e-3,
      //       2.7344809236135841598530720476098e-1,
      //       5.9764454459057452955269360703495e-3,
      //       5.2410938325006369651660186410547e-1,
      //       1.0057301050542603016660734438860e-2,
      //       9.0895615934244191891990428663041e-2,
      //       8.8452523398346606829848673809009e-1,
      //       3.8799010163202809566473695352386e-2,
      //       5.1227366986894646065998675430995e-1,
      //       3.9602679189978358269945323907258e-1,
      //       2.6174937160968349701196686726870e-2,
      //       7.6203914975198293014324875248846e-1,
      //       1.5901928319293625479920526665605e-1,
      //       3.9575280522419545364193216395460e-2,
      //       1.9753738413334949871287005642309e-1,
      //       2.5612186241563726204876368462414e-1,
      //       2.3188161972795448597914015255621e-2,
      //       6.8579014955455508229105658625928e-1,
      //       4.5690890104473341999976206420424e-2
      //     };

      // A degree=10 rule with 27 QPs
      // if (r.d==10 && r.nc==0 && r.nv==1 && r.ne==3 && r.ng==5)
      //   x =
      //     {
      //       7.0349169304985602109942089161121e-4,
      //       4.2428107036736195296584198917347e-3,
      //       8.9246341808720256261073890380435e-1,
      //       7.0016114161227155365721513157167e-3,
      //       3.0494756762697497590213637925978e-1,
      //       6.5757487042465969693037379103481e-3,
      //       6.0349372376167113477272117355728e-1,
      //       4.5367047076144017880509249966141e-2,
      //       2.7702918680309928564466502274207e-1,
      //       4.8213298384380707906450631732514e-1,
      //       3.3434147221659577726825933401689e-2,
      //       4.5537573281654840359879597498098e-1,
      //       8.1628735765941535287742569977895e-2,
      //       2.7372947537570407982425816418687e-2,
      //       6.8321697249249291081054474783040e-2,
      //       2.1461041150208394379962799554115e-1,
      //       2.8008115730342111834738616377667e-2,
      //       9.5199220006813800448433295908254e-2,
      //       6.7380547129004477728227558622451e-1,
      //       1.3960746583857763185533320493072e-2,
      //       8.6097589069840721814036297772080e-2,
      //       3.8205053391734468157969505461413e-2
      //     };

      // std::cout << "x=" << std::endl;
      // print(x);

      // Will be set to the minimum objective value, upon return.
      double minf;
      nlopt_result res = nlopt_optimize(opt, x.data(), &minf);

      // Call the optimization routine
      if (res < 0) {
        std::cout << "nlopt failed!\n";
      }
      else
        {
          std::cout << "nlopt converged with reason: " << nlopt_result_to_string(res) << std::endl;
          std::cout << "found minimum " << minf << " at " << std::endl;
          print(x);

          // Set solver_data.u = x
          solver_data.u.resize(dim);
          for (unsigned int i=0; i<dim; ++i)
            solver_data.u[i] = x[i];

          // Call our Newton solver, using the nlopt solution as the initial guess.
          // One really frustrating thing that I found in the d==7 case was that,
          // even when started very close to the true solution by using an initial
          // guess converged by nlopt, Newton could still diverge.
          bool converged = false;
          converged = newton_min(solver_data);
          converged = newton(solver_data);
          // converged = gradient_descent(solver_data);
          // converged = nlcg(solver_data);

          if (!converged)
            std::cout << "Solution is *not* converged." << std::endl;
          else
            {
              std::cout << "Solution is converged!" << std::endl;
              print(solver_data.u);
            }
        }

      // Debugging: exit after one iteration.
      // break;
    } // end while(true)

  // Clean up
  nlopt_destroy(opt);

  return 0;
}


double myfunc(unsigned n, const double * x, double * grad, void * my_func_data)
{
  auto solver_data = static_cast<SolverData *>(my_func_data);

  // Get references to residual, Jacobian, and solution vector stored on the
  // SolverData object.
  Matrix<mpfr_class> & jac = solver_data->jac;
  std::vector<mpfr_class> & u = solver_data->u;
  std::vector<mpfr_class> & r = solver_data->r;

  // Copy data. It would be better if we didn't have to do this. One idea
  // would be to make the SolverData class templated on numerical type...
  u.assign(x, x+n);
  solver_data->ro3.residual_and_jacobian(&r, grad ? &jac : nullptr, u);

  // The gradient is grad_f = J^T * r
  if (grad)
    {
      std::vector<mpfr_class> grad_f =
        jac.matvec_transpose(r);
      for (unsigned int i=0; i<n; ++i)
        grad[i] = grad_f[i].get_d();
    }

  // The objective function value is 0.5 * dot(r,r)
  return 0.5 * dot(r,r).get_d();
}



double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
  // x-index of orbit to which are applying the constraint.
  auto kptr = static_cast<unsigned int *>(data);
  auto k = *kptr;

  if (grad)
    {
      grad[k] = 1.;
      grad[k+1] = 1.;
    }
  return x[k] + x[k+1] - 1;
}


void usage()
{
  std::cout << "Usage: ./drivers/nlopt_test -d 10 -c 0 -v 1 -e 3 - g 5" << std::endl;
  std::cout << "-d = polynomial degree of exactness" << std::endl;
  std::cout << "-c = number of centroid orbits" << std::endl;
  std::cout << "-v = number of vertex orbits" << std::endl;
  std::cout << "-e = number of edge orbits" << std::endl;
  std::cout << "-g = number of general orbits" << std::endl;
}


std::string nlopt_result_to_string(nlopt_result res)
{
  switch (res)
    {
    case 1: // generic success code
      return std::string("NLOPT_SUCCESS");
    case 2:
      return std::string("NLOPT_STOPVAL_REACHED");
    case 3:
      return std::string("NLOPT_FTOL_REACHED");
    case 4:
      return std::string("NLOPT_XTOL_REACHED");
    case 5:
      return std::string("NLOPT_MAXEVAL_REACHED");
    case 6:
      return std::string("NLOPT_MAXTIME_REACHED");
    default:
      return std::string("Unrecognized nlopt_result " + std::to_string(res));
    }
}

#endif // HAVE_NLOPT
