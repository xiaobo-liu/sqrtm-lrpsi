function [X,iter] = sqrtm_in(A, scal, tol, maxit, mode)
% SQRTM_IN    Matrix square root iteration.
%    [X,ITER] = SQRTM_IN(A) is the principal square root of the m-by-m matrix A
%    A after ITER steps of the IN iteration
%
%    F_i = (E_i + X_i / 2) / mu_i - mu_i X_i / 2,
%    X_i = mu_i X_i + F_i,                             X_0 = A,          (*)
%    E_i = - (F_i X_i^{-1} F_i) / 2,                   E_0 = I,
%
%    [...] = SQRTM_IN(A,SCAL) uses the scaled version of the iteration if SCAL
%    is TRUE and the unscaled iteration otherwise. The default value of this
%    parameter is FALSE.
%
%    [...] = SQRTM_IN(A,SCAL,TOL) stops the iteration when the difference
%    between two successive iterates in the 1-norm is below TOL. The default
%    value of this parameter is EPS(class(A)).
%
%    [...] = SQRTM_IN(A,SCAL,TOL,MAXIT) performs at most MAXIT iterations, but
%    terminates earlier if the stopping criterion based on the paramter TOL is
%    met. The default value of this parameter is 50.
%
%   [...] = SQRTM_IN(A,SCAL,TOL,MAXIT,MODE) selects the code to be used to
%   evaluate the expression (F_i X_i^{-1} F_i) appearing in (*). The possible
%   choices are:
%     * MODE = 0 => E = F * inv(X) * F
%     * MODE = 1 => E = F * (X \ F)
%     * MODE = 2 => E = (F /X) * F)

  [m,n] = size(A);
  assert(m ==  n);

  if nargin < 2
    scal = false;
  end
  if nargin < 3
    tol = eps(class(A));
  else
    assert(tol > 0 && tol < 1);
  end
  if nargin < 4
    maxit = 50;
  else
    assert(maxit > 0 && round(maxit) == maxit);
  end
  if nargin < 5
    mode = 1;
  else
    assert(mode >= 0 && mode <= 2 && round(mode) == mode)
  end

  X = A;
  E = (eye(m) - A) / 2;

  for iter = 1:maxit
    if scal
      mu = abs(det(X)^2/det(A))^(-1/(2*m));
      if isnan(mu) || isinf(mu) || mu == 0
        mu = 1;
        scal = 0;
      end
    else
      mu = 1;
    end
    Etilde = (E+X/2)/mu - mu*X/2;
    Xnew = mu*X + Etilde;
    switch mode
      case 0
        E = -(Etilde * inv(Xnew) *Etilde) / 2;
      case 1
        E = -(Etilde * (Xnew \ Etilde)) / 2;
      case 2
        E = -((Etilde / Xnew) * Etilde) / 2;
    end
    if norm(X - Xnew, 1) / norm(X, 1) < tol
      break;
    end
    X = Xnew;
  end
end