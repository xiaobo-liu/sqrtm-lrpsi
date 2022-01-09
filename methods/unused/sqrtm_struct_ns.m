function [Y, iter] = sqrtm_struct_ns(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using a structured version of the
% Newton-Schulz iteration. The function does not check whether the sufficient
% condition for the convergence of the iteration, that is, NORM(I - A, P) < 1,
% for P = '1', '2', or 'Inf', is satisfied. SCAL specifies whether scaling is
% used (SCAL = 1) or not (SCAL = 0, default).

% The cost of the first iteration could be reduced, in principle, by taking into
% account the fact that D = 0. We don't do that for the moment.

% Y = gamma*eye(n) + U*C*V'
% Z = delta*eye(n) + U*D*V'
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

  % First few iterations can be done implicitly as no computation is required.
  % *** Iteraktion 0
  % gamma = epsilon;
  % C = eye(k);
  % delta = 1;
  % D = zeros(k);
  % *** Iteration 1 % implicit, without scaling
  delta = (3-epsilon)/2;
  gamma = epsilon*delta;
  Ik = eye(k);
  C = (3/2-epsilon) * Ik - 1/2 * M;
  D = -1/2 * Ik;

  % A = epsilon*eye(n) + U*V';
  % Y = gamma*eye(n) + U*C*V';
  % Z = delta*eye(n) + U*D*V';
  % norm(Y - A)
  % norm(Z - eye(n))

  if scal
      detA = abs(sqrt( det(epsilon*Ik + M) ))^(1/n);
      pre_mu = 1;
  end


  for iter = 1:maxit

    % Ynew = Y*(3*eye(n) -  Z*Y)/2;
    % Z = (3*eye(n) - Z*Y)*Z/2;
    % Y = Ynew;

    P = D * M;
    Q = P * C;
    S = gamma*D + delta*C + gamma*delta*Q;

    % Determinantal scaling.
    SM = S*M;
    if scal
        mu = (sqrt(epsilon)/gamma)^(1-k/n) * detA / abs(det(gamma*Ik + M*C))^(1/n);
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

    gamma_new = mu*gamma*(3-mu^2*gamma*delta)/2;
    C_new = (mu/2) * ((3-mu^2*gamma*delta)*C - mu^2*gamma*S - mu^2*C*M*S);

    delta_new = mu*delta*(3-mu^2*gamma*delta)/2;
    D_new = (mu/2) * ((3-mu^2*gamma*delta)*D - mu^2*delta*S - mu^2*SM*D);

    if norm(C-C_new, 1) / norm(C, 1) < tol
      break
    end

    gamma = gamma_new;
    C = C_new;
    delta = delta_new;
    D = D_new;

    % disp('******')
    % Yt = gamma*eye(n) + U*C*V';
    % Zt = delta*eye(n) + U*D*V';
    % norm(Y - Yt)
    % norm(Z - Zt)
    % X = alpha*eye(n) + U*C*V';
    % norm(X^2-A) / norm(A)

  end

  % Y  = alpha*eye(n) + U*C*V';
  Y = U*C_new*V';
  Y(1:n+1:n^2) = Y(1:n+1:n^2) + gamma_new;

end