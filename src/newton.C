// Local includes
#include "solver_data.h"
#include "newton.h"
#include "vect.h"
#include "ro3.h"

bool newton(SolverData & solver_data)
{
  // Newton iteration parameters.
  mpfr_class tol = solver_data.tol;
  mpfr_class divtol = solver_data.divtol;
  mpfr_class alphamin = solver_data.alphamin;
  bool do_backtracking = solver_data.do_backtracking;
  // If do_backtracking is true, do we also require the backtracking
  // to produce a residual reduction before accepting the step?
  bool residual_reduction_required = solver_data.residual_reduction_required;
  unsigned iter = 0;
  unsigned int maxits = solver_data.maxits;
  bool converged = false;

  Ro3 & ro3 = solver_data.ro3;

  std::vector<mpfr_class> & u = solver_data.u;
  std::vector<mpfr_class> & r = solver_data.r;
  Matrix<mpfr_class> & jac = solver_data.jac;

  unsigned int n = u.size();

  // Storage for residual, Newton update, and Jacobian
  std::vector<mpfr_class> du(n), trial_u(n);

  // Store previous residual norm. Used for simplified line searching technique.
  // mpfr_class old_residual_norm = 0.;

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual (vector) at the current guess.
      ro3.residual_and_jacobian(&r, nullptr, u);

      // Check the norm of the residual vector to see if we are done.
      mpfr_class residual_norm = norm(r);

      if (solver_data.verbose)
        std::cout << "newton: Iteration " << iter
                  << ", residual_norm=" << residual_norm
                  << std::endl;

      if (residual_norm < tol)
        {
          converged = true;
          break;
        }

      // If the residual is too large, just give up.
      if (residual_norm > divtol)
        break;

      // Compute Jacobian (only) at the current guess.
      ro3.residual_and_jacobian(nullptr, &jac, u);

      // Compute update: du = -jac^{-1} * r
      try
        {
          jac.lu_solve(du, r);
        }
      catch (const std::exception & e)
        {
          // The Jacobian is singular, so we didn't converge. Break out of
          // for-loop with converged==false.
          std::cout << "newton(): " << e.what() << std::endl;
          break;
        }

      // Compute next iterate, u -= alpha*du using simplified
      // backtracking, alpha_min < alpha <= 1.
      // This is only theoretically helpful at the moment, so far I have not encountered
      // an actual case where a failing case was able to converge by using backtracking.
      if (do_backtracking)
        {
          mpfr_class alpha = mpfr_class(1);
          bool linesearch_converged = false;
          // 2^(-50) ~ 9e-16
          for (unsigned int back=0; back<25; ++back)
            {
              // std::cout << "Trying step with alpha=" << alpha << std::endl;
              trial_u = u - alpha * du;
              ro3.residual_and_jacobian(&r, nullptr, trial_u);
              mpfr_class trial_residual = norm(r);

              // if (solver_data.verbose)
              //   std::cout << "alpha = " << alpha
              //             << ", trial_residual = " << trial_residual << std::endl;

              // Was the residual reduced? (Or do we actually care?)
              // If we don't care, just set this flag to true.
              bool residual_reduced =
                residual_reduction_required ? (trial_residual < residual_norm) : true;

              // Check feasibility of the proposed solution.
              bool feasible = ro3.check_feasibility(trial_u);

              // if (solver_data.verbose && !residual_reduced)
              //   std::cout << "newton: backtracking step "
              //             << back
              //             << ", Residual was not reduced!"
              //             << std::endl;

              // if (solver_data.verbose && !feasible)
              //   std::cout << "newton: backtracking step "
              //             << back
              //             << ", trial solution was not feasible!"
              //             << std::endl;


              if (residual_reduced && feasible)
                {
                  // Accept this trial solution
                  u = trial_u;
                  linesearch_converged = true;

                  // if (solver_data.verbose)
                  //   std::cout << "Accept step with alpha = " << alpha << std::endl;

                  break;
                }
              else
                alpha /= 2;
            }

          // If the linesearch failed, we can quit or just keep going.
          if (!linesearch_converged)
            {
              // 1.) Keep going
              // subtract_scaled(u, alpha, du);

              // 2.) Quit
              std::cout << "Linesearch failed to find acceptable step." << std::endl;
              break;
            }
        }
      else
        {
          // Accept full step
          u -= du;
        }
    } // end while (true)

  // if (!converged)
  //   std::cout << "Newton iteration diverged, backtracking failed, or max iterations exceeded."
  //             << std::endl;

  return converged;
}



bool newton_min(SolverData & solver_data)
{
  // Newton minimization iteration parameters.
  mpfr_class tol = solver_data.tol;
  mpfr_class divtol = solver_data.divtol;
  mpfr_class alphamin = solver_data.alphamin;
  unsigned int maxits = solver_data.maxits;
  bool do_backtracking = solver_data.do_backtracking;
  bool residual_reduction_required = solver_data.residual_reduction_required;
  unsigned iter = 0;
  bool converged = false;

  Ro3 & ro3 = solver_data.ro3;

  std::vector<mpfr_class> & u = solver_data.u;
  std::vector<mpfr_class> & r = solver_data.r;
  Matrix<mpfr_class> & jac = solver_data.jac;

  // Problem size
  unsigned int n = u.size();

  // Storage for residual, Newton update, and Jacobian
  std::vector<mpfr_class> du(n), grad_f(n), trial_u(n);
  Matrix<mpfr_class> djac(n,n), hess(n,n);

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual and Jacobian at the current guess.
      ro3.residual_and_jacobian(&r, &jac, u);

      // Compute the size of the scalar function "f"
      // which we are trying to minimize, _not_ dot(r,r)^0.5.
      mpfr_class residual = 0.5 * dot(r,r);

      if (solver_data.verbose)
        std::cout << "newton_min: Iteration " << iter
                  << ", residual = " << residual
                  << std::endl;

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

      // if (solver_data.verbose)
      //   std::cout << "|grad_f| = " << norm(grad_f) << std::endl;

      // Compute the Hessian. We will approximate the true Hessian by
      // finite differencing the Jacobian.
      hess.clear();
      hess.resize(n,n);

      // One part of the Hessian is jac^T * jac. I think this is guaranteed to be SPD?!
      for (unsigned int i=0; i<n; ++i)
        for (unsigned int j=0; j<n; ++j)
          for (unsigned int k=0; k<n; ++k)
            hess(i,j) += jac(k,i) * jac(k,j);

      // Debugging: check symmetry of hessian
      // if (solver_data.verbose)
      //   std::cout << "J^T * J asymmetry=" << hess.asymmetry() << std::endl;

      // Finite differencing parameter. I'm not sure I actually want
      // to include this. J^T*J is guaranteed to be SPD, if we add
      // this bit to it, that is no longer guaranteed. Also it's very
      // expensive to do the finite differencing, and finally, the
      // term we are adding is small when the residual is already
      // small.
      //
      // Note: I confirmed that adding this contribution to the Hessian maintains
      // its symmetry, but apparently makes it not SPD? Whereas J^T * J is always
      // SPD... One idea is that eps should probably depend on how close we are
      // to the minimum, one way to do that is to have it be proportional to |grad_f|.
//      mpfr_class eps = mpfr_class(1.e-8) * norm(grad_f);
//      for (unsigned int j=0; j<n; ++j)
//        {
//          std::vector<mpfr_class> u_eps = u;
//          u_eps[j] += eps;
//
//          // Compute Jacobian at u + eps
//          ro3.residual_and_jacobian(nullptr, &djac, u_eps);
//
//          // Debugging
//          // std::cout << "j=" << j << std::endl;
//          // djac.print();
//
//          for (unsigned int k=0; k<n; ++k)
//            for (unsigned int i=0; i<n; ++i)
//              hess(i,j) += r[k] * (djac(k,i) - jac(k,i)) / eps;
//        }

//      if (solver_data.verbose)
//        std::cout << "Full Hessian asymmetry=" << hess.asymmetry() << std::endl;

      // Debugging: print Hessian.
      // std::cout << "hess=" << std::endl;
      // hess.print();

      // Solve Hessian system for update, -du.
      while (true)
        {
          try
            {
              auto hess_copy = hess;
              // Note: it seems to be possible to get different
              // solutions depending on if you solve with Cholesky or
              // LU...
              hess_copy.cholesky_solve(du, grad_f);
              // hess_copy.lu_solve(du, grad_f);
              break;
            }
          catch (const std::exception & e)
            {
              // Matrix is singular or not SPD. Let's add beta*I
              // std::cout << "newton_min(): " << e.what() << std::endl;

              // We can add some a factor which is proportional to the
              // diagonal existing diagonal but that didn't work very well.
              mpfr_class beta(1);
              for (unsigned int i=0; i<n; ++i)
                hess(i,i) += beta; // * abs(hess(i,i));
            }
        }

      // Debugging
      // if (solver_data.verbose)
      //   {
      //     std::cout << "du=" << std::endl;
      //     print(du);
      //   }

      // Line search
      mpfr_class alpha = mpfr_class(1);
      bool linesearch_converged = false;

      if (do_backtracking)
        {
          for (unsigned int back=0; back<50; ++back)
            {
              // std::cout << "Trying step with alpha=" << alpha << std::endl;
              trial_u = u - alpha * du;

              // Compute residual at trial solution.
              ro3.residual_and_jacobian(&r, nullptr, trial_u);
              mpfr_class trial_residual = 0.5 * dot(r,r);

              // Check for feasibility of the solution
              bool feasible = ro3.check_feasibility(trial_u);
              bool residual_reduced =
                residual_reduction_required ? (trial_residual < residual) : true;

              // if (solver_data.verbose)
              //   std::cout << "feasible = " << feasible << std::endl;

              if (!feasible || !residual_reduced)
                {
                  alpha /= 2;

                  // Debugging:
                  // if (solver_data.verbose)
                  //   std::cout << "Trying step with alpha=" << alpha << std::endl;
                }
              else
                {
                  // Accept this trial solution
                  linesearch_converged = true;
                  break;
                }
            }
        }

      if (!linesearch_converged)
        {
          std::cout << "Linesearch failed to find acceptable step." << std::endl;
          break;
        }

      subtract_scaled(u, alpha, du);

      // Debugging: print du
      // std::cout << "du=" << std::endl;
      // print(du);
    } // end while (true)

  return converged;
}
