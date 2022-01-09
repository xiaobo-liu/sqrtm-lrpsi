function Y = sqrtm_unstruct_ns(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using the Newton-Schulz iteration.
% The function does not check whether the sufficient condition for the
% convergence of the iteration, that is, NORM(I - A, P) < 1, for P = '1', '2',
% or 'Inf', is satisfied.

  if nargin < 5
    tol = eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end

  m = size(U, 1);
  A = epsilon*eye(m) + U*V';
  m = size(A, 1);
  Y = A;
  Z = eye(m);

  for i = 1:maxit
    if scal
      mu = abs(det(X)^2/det(A))^(-1/(2*m));
    else
      mu = 1;
    end
    Ynew = mu*Y*(3*eye(m) -  mu^2*Z*Y)/2;
    Z = mu*(3*eye(m) - mu^2*Z*Y)*Z/2;
    Y = Ynew;

    if norm(Y-Ynew, 1) / norm(Ynew, 1) < tol
      Y = Ynew;
      break
    end
    Y = Ynew;
  end

end