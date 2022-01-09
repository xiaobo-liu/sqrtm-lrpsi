function [Y, iter] = sqrtm_struct_db_prod(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using a structured version of
% the Denman-Beavers iteration in product form. SCAL specifies whether
% scaling is used (SCAL = 1) or not (SCAL = 0, default).

  [n,k] = size(U);

  if nargin < 4
    scal = false;
  end
  if nargin < 5
    tol = eps(class(U));
  end
  if nargin < 6
    maxit = 50;
  end


  M = V' * U;

  nu = epsilon;
  gamma = epsilon;
  Ik = eye(k);
  N = Ik;
  C = Ik;


  % A = epsilon*eye(n) + U*V';
  % W = nu*eye(n) + U*N*V';
  % Y = gamma*eye(n) + U*C*V';

  if scal, pre_mu = 1; end

  for iter = 1:maxit

    MN = M*N;
    % MN(1:k+1:k^2) = MN(1:k+1:k^2)+ nu;
    MN = MN + nu*Ik;

    % Determinantal scaling.
    if scal
      mu = sqrt( nu^(1-k/n) * abs(det(MN))^(1/n));
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

    S = (N/MN) / (mu^2*nu);

    nu_new = (1 + ( mu^2*nu+1/(mu^2*nu) )/2 )/2;
    N_new = (mu^2*N - S)/4;

    gamma_new = ( mu*gamma + gamma/(mu*nu)) /2;
    C_new = ((mu + 1/(mu*nu))*C - mu*gamma*S - mu*C*M*S )/ 2;

    if norm(C-C_new, 1) / norm(C, 1) < tol
      break
    end

    nu = nu_new;
    N = N_new;
    gamma = gamma_new;
    C = C_new;
  end

  % Y  = gamma*eye(n) + U*C*V';
  Y = U*C_new*V';
  Y(1:n+1:n^2) = Y(1:n+1:n^2) + gamma_new;

end
