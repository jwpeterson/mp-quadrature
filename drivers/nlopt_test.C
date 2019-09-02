// C++ includes
#include <iostream>

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

// nlopt includes
#include <nlopt.h>

// C/C++ includes
#include <math.h>
#include <stdlib.h> // random()

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
  auto solver_data = static_cast<SolverData *>(my_func_data);

  // Copy data. It would be better if we didn't have to do this. One idea
  // would be to make the SolverData class templated on numerical type...
  std::vector<mpfr_class> r;
  Matrix<mpfr_class> jac;
  std::vector<mpfr_class> u(x, x+n);
  solver_data->residual_and_jacobian(&r, grad ? &jac : nullptr, u);

  // The gradient is grad_f = J^T * r
  if (grad)
    {
      std::vector<mpfr_class> grad_f = jac.matvec_transpose(r);
      for (unsigned int i=0; i<n; ++i)
        grad[i] = grad_f[i].get_d();
    }

  // The objective function value is 0.5 * dot(r,r)
  return 0.5 * dot(r,r).get_d();
}


// Constraints are of the form fc(x) <= 0, where fc is required to have
// the same parameters as the objective function. In our case, there is
// one inequality constraint associated with each ro3 orbit of the form:
// 1 - x - y >= 0,
// which we multiply by -1 to write as:
// x + y - 1 <= 0
// So in this case, the incoming "data" pointer just needs to tell us
// the index of the x dof for the orbit being constrained.
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

// Convert a numerical nlopt_result successful return code to a human readable string.
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

int main()
{
  std::cout.precision(16);
  std::cout.setf(std::ios_base::scientific);

  // Use the current time since epoch as a seed.
  time_t seed = time(nullptr);
  // long seed = 1567251384;
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  unsigned int d=10;

  SolverData solver_data;
  solver_data.verbose = true;
  solver_data.maxits = 50;
  solver_data.residual_and_jacobian = ResidualAndJacobian(d);
  solver_data.check_feasibility = CheckFeasibility(d);
  solver_data.do_backtracking = true;
  solver_data.residual_reduction_required = true;

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

  // The problem dimension depends only on "d".
  unsigned int dim = (d*d + 3*d + 6)/6;

  // There will be no centroid if dim is a multiple of 3.
  unsigned int n_centroid = dim % 3;

  // The number of Ro3 orbits
  unsigned int n_ro3 = dim / 3;

  if (n_centroid == 2)
    {
      std::cout << "Error: rules with 2 centroid points don't make sense." << std::endl;
      exit(1);
    }

  // Set minimization algorithm and dimensionality
  nlopt_opt opt = nlopt_create(alg, dim);

  if (alg == NLOPT_AUGLAG)
    {
      // The objective function, bounds, and nonlinear-constraint
      // parameters of local_opt are ignored.
      nlopt_opt local_opt = nlopt_create(NLOPT_LD_SLSQP, dim);
      nlopt_set_xtol_rel(local_opt, 1.e-8);
      nlopt_set_ftol_rel(local_opt, 1.e-16);
      nlopt_set_local_optimizer(opt, local_opt);
    }

  // lower bounds: all parameters must be positive, so we just use the
  // default-constructed values of the std vector for lb. All weights
  // for 3-point orbits must be less than 1/6, which is (1/3)*(ref. element area).
  // Everything else must be less than 1.
  std::vector<double> lb(dim);
  std::vector<double> ub;

  double sixth = 1./6;
  if (n_centroid)
    ub.push_back(sixth);

  for (unsigned int i=0; i<n_ro3; ++i)
    {
      ub.push_back(sixth);
      ub.push_back(1);
      ub.push_back(1);
    }

  nlopt_set_lower_bounds(opt, lb.data());
  nlopt_set_upper_bounds(opt, ub.data());
  nlopt_set_min_objective(opt, myfunc, static_cast<void *>(&solver_data));

  // Add inequality constraint. We need call the "add" function once
  // per inequality constraint, and each time pass it the index of the
  // x DOF of the orbit which is to be constrained.
  std::vector<unsigned int> data;
  for (unsigned int i=0; i<n_ro3; ++i)
    data.push_back(3*i + 1 + n_centroid);

  // print(data);

  for (unsigned int c=0; c<data.size(); ++c)
    nlopt_add_inequality_constraint(opt, myconstraint, &data[c], 1e-8);

  // Set tolerance(s)
  nlopt_set_xtol_rel(opt, 1.e-16);
  nlopt_set_ftol_rel(opt, 1.e-16);

  // You can use this parameter to control how long some of the
  // gradient free optimization algorithms run, otherwise it seems
  // that they can just run forever.
  // nlopt_set_maxeval(opt, 9999);

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

  // Set initial guess.
  // std::vector<double> x =
  //   {
  //     2.65e-02, // w1
  //     6.23e-02, // x1
  //     6.75e-02, // y1
  //
  //     4.38e-02, // w2
  //     5.52e-02, // x2
  //     3.21e-01, // y2
  //
  //     2.87e-02, // w3
  //     3.43e-02, // x3
  //     6.60e-01, // y3
  //
  //     6.74e-02, // w4
  //     5.15e-01, // x4
  //     2.77e-01, // y4
  //   };

  while (true)
    {
      std::vector<double> x;
      if (n_centroid)
        x.push_back(sixth * double(random())/RAND_MAX);

      for (unsigned int i=0; i<n_ro3; ++i)
        {
          x.push_back(sixth * double(random())/RAND_MAX);
          x.push_back(double(random())/RAND_MAX);
          x.push_back(double(random())/RAND_MAX);
        }

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
          bool converged = newton(solver_data);

          if (!converged)
            std::cout << "Solution is *not* converged." << std::endl;
          else
            {
              std::cout << "Solution is converged!" << std::endl;
              print(solver_data.u);
            }
        }
    } // end while(true)

  // Clean up
  nlopt_destroy(opt);

  return 0;
}

#endif // HAVE_NLOPT
