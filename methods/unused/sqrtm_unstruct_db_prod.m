function X = sqrtm_unstruct_db_prod(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using the DB iteration in 
% product form.

  if nargin < 5
    tol = eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end

  m = size(U,1);
  X = epsilon*(eye(m) + U*V');
  M = X;

  for i=1:maxit
      
      if scal
          mu = abs(det(M))^(-1/(2*m));
          X = mu*X; M = mu^2*M;
      end
      invM = inv(M);
      
      Xnew = X*(eye(m) + invM)/2;
      M = 0.5*(eye(m) + (M + invM)/2);
    
      if norm(Xnew-X, 1) / norm(Xnew, 1) < tol
          X = Xnew;
          break;
      end
      X = Xnew;
  end

end