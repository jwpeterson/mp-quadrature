// mp-quadrature includes
#include "solver_data.h"
#include "newton.h"

// nlopt includes
#include <nlopt.h>

// C/C++ includes
#include <math.h>
#include <iostream>
#include <stdlib.h> // random

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

  // On the first random seed attempted with AUGLAG + NLOPT_LD_SLSQP, I got an answer!
  // And in fact it converges for almost every random seed I give it!
  // seed=1567253290
  // found minimum 2.1608544652295794e-30 at
  // 6.7493187010515404e-02
  // 2.0644149863649774e-01
  // 2.7771616701986818e-01
  // 2.6517028155714859e-02
  // 6.7517867067935791e-02
  // 6.2382265096082017e-02
  // 2.8775042779261711e-02
  // 3.0472650088237996e-01
  // 6.6094919618974801e-01
  // 4.3881408721174679e-02
  // 5.5225456672372955e-02
  // 6.2327204951356407e-01

  // NLOPT_LD_TNEWTON_PRECOND_RESTART
  // 2.6415309844487228e-02
  // 6.2048514582985367e-02
  // 6.7528655243528232e-02
  // 4.4133843200339001e-02
  // 5.6098712400163488e-02
  // 3.2069483854069214e-01
  // 2.9725408563411330e-02
  // 3.5408617008180480e-02
  // 6.5952241159895886e-01
  // 6.6392105258589698e-02
  // 5.1441538853983093e-01
  // 2.7635340137105119e-01

  // NLOPT_LD_SLSQP - this didn't really converge when I started right
  // at the known solution, but it did converge from a random initial
  // condition created from this seed...
  // seed=1567251384
  // found minimum 1.2735061228156008e-28 at
  // 2.8775042826050961e-02
  // 6.6094919599171909e-01
  // 3.0472650096118875e-01
  // 6.7493186999977570e-02
  // 5.1584233436403581e-01
  // 2.7771616681035760e-01
  // 4.3881408676599953e-02
  // 5.5225456584071402e-02
  // 3.2150249387976276e-01
  // 2.6517028164038121e-02
  // 6.2382265127586128e-02
  // 6.7517867051875763e-02

  // NLOPT_LD_MMA, xtol_rel=1.e-8, objective=1.6210802277134142e-15
  // 2.6569441070361466e-02
  // 6.2570257228477255e-02
  // 6.7453629074729418e-02
  // 4.3865680847562004e-02
  // 5.5424307869906826e-02
  // 3.2097324212345357e-01
  // 2.8767507727108708e-02
  // 3.4109356620012313e-02
  // 6.6010124639009604e-01
  // 6.7464037019301912e-02
  // 5.1580657331704916e-01
  // 2.7767321619059238e-01

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
  }

  // Clean up
  nlopt_destroy(opt);

  return 0;
}

// For degree=8, the AUGLAG+SLSQP approach finds the following results.
// They don't appear to be a "true" solutions, because they don't obtain
// a minium of ~1.e-28 like we see in other cases.
//
// seed=1567254545
// found minimum 4.1733578212990206e-14 at
// 5.1288152148630424e-02
// 3.2798183971138672e-01
// 6.0839572272242337e-01
// 2.9386634760013963e-02
// 3.1056160616542028e-01
// 3.5298840762582223e-01
// 6.0455344995865107e-02
// 9.3328463220219787e-02
// 6.2691829094939000e-01
// 7.2502902552967763e-04
// 8.3207209900354839e-01
// 3.2554053036089475e-01
// 2.4811506271759692e-02
// 8.7487291135983447e-01
// 5.6679536910645145e-02
//
// seed=1567254758
// found minimum 3.4895214406050540e-15 at
// 4.8819978580307167e-03
// 9.4930065786400064e-01
// 5.9772629477515715e-02
// 4.4482224069351164e-02
// 6.6931660863610049e-02
// 6.6876847799853989e-01
// 7.2437309288881452e-02
// 3.3190135224045408e-01
// 4.9413289458984316e-01
// 2.8653848858907319e-02
// 1.3286867079073769e-01
// 8.1593729193429132e-01
// 1.6211286584277510e-02
// 5.7678833060036661e-01
// 7.2785418952812011e-03
//
// seed=1567256705
// found minimum 1.0540286142691775e-14 at
// 4.1145098816673602e-02
// 1.8982794358987251e-01
// 1.1851008643947880e-01
// 1.6515284424203813e-02
// 5.5280049840639653e-02
// 4.6016215661939773e-02
// 6.9846820440213231e-02
// 3.7533297695222645e-01
// 4.6193896193529088e-01
// 2.0615070136278047e-02
// 1.6929612202094126e-02
// 2.7355278183716908e-01
// 1.8544393143362072e-02
// 6.2390091356454391e-01
// 3.6140865438068071e-01
//
// seed=1567256852
// found minimum 6.0374269209491919e-15 at
// 1.2962322841122963e-02
// 0.0000000000000000e+00
// 5.1309571799085185e-01
// 1.2399095040475129e-02
// 5.5929055887210061e-02
// 3.1749396125371872e-02
// 4.4970568329451927e-02
// 6.8528679877422272e-01
// 2.3379935293889950e-01
// 7.1983762206953378e-02
// 1.6726108419373828e-01
// 3.5313980074743351e-01
// 2.4350918287324667e-02
// 2.0657251688380829e-01
// 7.5787844664414283e-01
//
// d=10
// seed=1567284954
// nlopt converged with reason: NLOPT_FTOL_REACHED
// found minimum 5.9417632242414094e-17 at
// 2.2718419969892670e-02
// 3.6550053800521981e-02
// 6.2929337143671449e-01
// 2.6739391552225406e-01
// 1.3487985098020129e-02
// 1.7221003064418045e-01
// 2.6305098047661093e-02
// 9.4610102916900338e-03
// 3.3141816740849067e-02
// 4.4334931361701801e-02
// 6.3007287326519931e-03
// 5.4662327958830503e-01
// 2.2901824033310489e-01
// 4.7199961255077416e-02
// 1.6502056063496212e-01
// 3.5635696686623852e-01
// 1.6853815288768500e-02
// 4.9883179267455252e-01
// 1.3405403244585761e-02
// 2.9240305547971579e-02
// 5.3756906280127539e-02
// 2.0598528906603203e-01
