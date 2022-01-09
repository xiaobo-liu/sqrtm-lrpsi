function [X,i] = sqrtm_unstruct_in(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using the IN iteration.

  if nargin < 4
    scal = false;
  end
  if nargin < 5
    tol = eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end
  m = size(U,1);

  A = epsilon*eye(m) + U*V';
  X = A;
  E = (eye(m) - A) / 2;

  for i = 1:maxit
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
    E = -(Etilde * (Xnew \ Etilde)) / 2;
    % E = -(Etilde * (inv(Xnew) * Etilde)) / 2;
    % condXnew = cond(Xnew)
    % norm(E - Enew, 1) / norm(E, 1)
    if norm(X - Xnew, 1) / norm(X, 1) < tol
      break;
    end
    X = Xnew;
  end
end