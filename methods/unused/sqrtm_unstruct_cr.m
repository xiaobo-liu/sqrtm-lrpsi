function Z = sqrtm_unstruct_cr(epsilon, U, V, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using the CR iteration.

  if nargin < 5
    tol = eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end
  m = size(U,1);

  A = epsilon*eye(m) + U*V';
  Y = eye(m) - A;
  Z = 2 * (A + eye(m));

  for i=1:maxit
    Y = -Y * (Z \ Y);
    Znew = Z + 2*Y;
    if norm(Znew-Z, 1) / norm(Znew, 1) < tol
      Z = Znew / 4;
      break;
    end
    Z = Znew;
  end
end