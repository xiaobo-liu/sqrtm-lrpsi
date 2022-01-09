function X = sqrtm_unstruct_db(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using the DB iteration.

  if nargin < 5
    tol = eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end

  m = size(U,1);
  X = epsilon*eye(m) + U*V';
  Y = eye(m);

  for i=1:maxit
    if scal
      mu = abs(det(X)^2/det(A))^(-1/(2*m));
    else
      mu = 1;
    end
    Xnew = (mu*X + inv(Y)/mu) / 2;
    Y = (mu*Y + inv(X)/mu) / 2;
    if norm(Xnew-X, 1) / norm(Xnew, 1) < tol
      X = Xnew;
      break;
    end
    X = Xnew;
  end

end