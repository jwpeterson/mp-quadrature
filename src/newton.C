// Local includes
#include "solver_data.h"
#include "newton.h"
#include "vect.h"

bool newton(SolverData & solver_data)
{
  // Newton iteration parameters.
  mpfr_class tol = solver_data.tol;
  mpfr_class divtol = solver_data.divtol;
  mpfr_class alphamin = solver_data.alphamin;
  const bool do_backtracking = false;
  unsigned iter = 0;
  unsigned int maxits = solver_data.maxits;
  bool converged = false;

  ResidualAndJacobianFunctionPtr
    residual_and_jacobian = solver_data.residual_and_jacobian;

  std::vector<mpfr_class> & u = solver_data.u;

  // Storage for residual, Newton update, and Jacobian
  std::vector<mpfr_class> r, du;
  Matrix<mpfr_class> jac;

  // Store previous residual norm. Used for simplified line searching technique.
  mpfr_class old_residual_norm = 0.;

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual (vector) at the current guess.
      (*residual_and_jacobian)(&r, nullptr, u);

      // Check the norm of the residual vector to see if we are done.
      mpfr_class residual_norm = norm(r);

      if (solver_data.verbose)
        std::cout << "Iteration " << iter << ", residual_norm=" << residual_norm << std::endl;

      if (residual_norm < tol)
        {
          converged = true;
          break;
        }

      // If the residual is too large, just give up.
      if (residual_norm > divtol)
        break;

      // Compute Jacobian (only) at the current guess.
      (*residual_and_jacobian)(nullptr, &jac, u);

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
      mpfr_class alpha(1);
      bool backtracking_converged = false;
      while (true)
        {
          if (alpha < alphamin)
            break;

          // This looks "nice" but is probably less efficient in general because
          // a temporary has to be created for the result of "alpha * du" and then
          // finally our custom operator-=() function is called. It's not possible
          // to define a non-member operator-=() that takes more than two args.
          // u -= alpha * du;

          // This looks less "cool" but is probably more efficient
          // because no temporary is created.
          subtract_scaled(u, alpha, du);

          if (do_backtracking)
            {
              (*residual_and_jacobian)(&r, nullptr, u);
              mpfr_class new_residual_norm = norm(r);

              if (solver_data.verbose)
                std::cout << "Trying Newton step with alpha = " << alpha
                          << ", residual = " << new_residual_norm << std::endl;

              if (new_residual_norm < residual_norm)
                {
                  backtracking_converged = true;
                  break;
                }

              // Don't waste time backtracking if we already diverged
              if (new_residual_norm > divtol)
                break;

              alpha /= mpfr_class(2);
            }
          else
            {
              // If we aren't doing backtracking, then backtracking is "done".
              backtracking_converged = true;
              break;
            }
        }

      if (!backtracking_converged)
        break;
    } // end while

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
  unsigned iter = 0;
  bool converged = false;

  ResidualAndJacobianFunctionPtr
    residual_and_jacobian = solver_data.residual_and_jacobian;

  std::vector<mpfr_class> & u = solver_data.u;

  // Problem size
  unsigned int n = u.size();

  // Storage for residual, Newton update, and Jacobian
  std::vector<mpfr_class> r(n), du(n), grad_f(n), trial_u(n);
  Matrix<mpfr_class> jac(n,n), djac(n,n), hess(n,n);

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual and Jacobian at the current guess.
      (*residual_and_jacobian)(&r, &jac, u);

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

      // Compute the Hessian. We will approximate the true Hessian by
      // finite differencing the Jacobian.
      hess.clear();
      hess.resize(n,n);

      // One part of the Hessian is jac^T * jac
      for (unsigned int i=0; i<n; ++i)
        for (unsigned int j=0; j<n; ++j)
          for (unsigned int k=0; k<n; ++k)
            hess(i,j) += jac(k,i) * jac(k,j);

      // Choosing beta=1 initially is equivalent to taking one
      // step of gradient descent. Note: this isn't part of the real
      // Hessian, but it may be necessary to ensure the Hessian remains
      // positive definite?
      // mpfr_class beta(1);
      // for (unsigned int i=0; i<n; ++i)
      //   hess(i,i) += beta;

      // Finite differencing parameter. I tried several different values
      // for this, but the results did not seem to be very sensitve to it?
      mpfr_class eps(1.e-8);

      for (unsigned int j=0; j<n; ++j)
        {
          std::vector<mpfr_class> u_eps = u;
          u_eps[j] += eps;

          // Compute Jacobian at u + eps
          (*residual_and_jacobian)(nullptr, &djac, u_eps);

          // Debugging
          // std::cout << "j=" << j << std::endl;
          // djac.print();

          for (unsigned int k=0; k<n; ++k)
            for (unsigned int i=0; i<n; ++i)
              hess(i,j) += r[k] * (djac(k,i) - jac(k,i)) / eps;
        }

      // Debugging: print Hessian.
      // std::cout << "hess=" << std::endl;
      // hess.print();

      // Solve Hessian system for update, -du.
      try
        {
          hess.lu_solve(du, grad_f);
        }
      catch (const std::exception & e)
        {
          // Matrix is singular, so we break out of the while loop
          // with converged==false.
          std::cout << "newton_min(): " << e.what() << std::endl;
          break;
        }

      // Debugging
      // print(du);

      // Line search
      mpfr_class alpha = mpfr_class(1);
      bool linesearch_converged = false;
      for (unsigned int back=0; back<10; ++back)
        {
          // std::cout << "Trying step with alpha=" << alpha << std::endl;
          trial_u = u - alpha * du;
          (*residual_and_jacobian)(&r, nullptr, trial_u);
          mpfr_class trial_residual = 0.5 * dot(r,r);

          // Debugging
          // std::cout << "alpha = " << alpha
          //           << ", trial_residual = " << trial_residual << std::endl;
          if (trial_residual < residual)
            {
              // Accept this trial solution
              u = trial_u;
              linesearch_converged = true;
              break;
            }
          else
            alpha /= 2;
        }

      // If the linesearch failed, stop iterating. Perhaps we could
      // instead take a gradient descent step?
      if (!linesearch_converged)
        {
          // It isn't actually terrible to allow a residual increase... since the
          // algorithm doesn't really have any "history" it's like it starts again
          // from a new initial guess and just keeps going.
          //
          // One thing I do see is when you encounter many linesearch failures in
          // a row, it usually means that no more progress is going to be made.
          // It might be good to detect this and stop iterating if this happens...

          // std::cout << "Linesearch failed to find a residual reduction." << std::endl;
          // break;

          // Accept the step even though it's not a residual reduction and keep going?
          subtract_scaled(u, alpha, du);
        }

      // Debugging: print du
      // std::cout << "du=" << std::endl;
      // print(du);
    } // end while (true)

  return converged;
}
