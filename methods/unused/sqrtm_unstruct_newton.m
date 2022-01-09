function X = sqrtm_unstruct_newton(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using the Newton's method.

  if nargin < 5
    tol = eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end
  m = size(U,1);
  A = epsilon*eye(m) + U*V';
  X = A;

  for i = 1:maxit
    if scal
      mu = abs(det(X)^2/det(A))^(-1/(2*m));
    else
      mu = 1;
    end
    Xnew = (mu*X + (mu*X)\A) / 2;
    if norm(Xnew-X, 1) / norm(Xnew, 1) < tol
      X = Xnew;
      break;
    end
  end

end