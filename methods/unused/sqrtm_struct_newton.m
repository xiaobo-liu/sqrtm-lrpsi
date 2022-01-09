function [X, iter] = sqrtm_struct_newton(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using a structured version of the
% Newton iteration. SCAL specifies whether scaling is used (SCAL = 1) or not
% (SCAL = 0, default).

% TODO: Can the cost of the first few iterations be reduced in this case?

  [n,k] = size(U);

  if nargin < 4
    scal = false;
  end
  if nargin < 5
    tol = 10*eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end


  M = V' * U;
  Ik = eye(k);

  beta = epsilon;
  B = Ik;

  % A = epsilon*eye(n) + U*V';
  % X = beta*eye(n) + U*B*V';
  % norm(X - A)

  if scal
      detA = abs(sqrt( det(epsilon*eye(k) + M) ))^(1/n);
      pre_mu = 1;
  end

  for iter = 1:maxit

    MB = M * B;
    % Determinantal scaling.
    if scal
        mu = (sqrt(epsilon)/beta)^(1-k/n) * detA / abs(det(beta*Ik + MB))^(1/n);
        if isnan(mu) || isinf(mu) || mu == 0
            mu = 1;
            scal = 0;
        end
        if abs((mu - pre_mu)/mu) < 0.01
            mu = 1;
            scal = 0;
        end
        pre_mu = mu;
    else
      mu = 1;
    end


    mubeta = mu*beta;
    beta_new = (mubeta + epsilon/mubeta) / 2;

    B_new = (mu*B + (Ik - B* ((beta*Ik+MB)\(epsilon*Ik+M)) )/mubeta ) / 2;

    if norm(B-B_new, 1) / norm(B, 1) < tol
      break
    end

    beta = beta_new;
    B = B_new;

    % disp('******')
    % Xt = beta*eye(n) + U*B*V';
    % norm(X - Xt)
    % norm(Xt^2-A) / norm(A)
    % norm(X^2-A) / norm(A)

    % keyboard

  end

  X = U*B_new*V';
  X(1:n+1:n^2) = X(1:n+1:n^2) + beta_new;

end
