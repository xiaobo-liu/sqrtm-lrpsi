function [X, iter] = sqrtm_struct_in(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using a structured version of
% the IN iteration. SCAL specifies whether scaling is used (SCAL = 1) or
% not (SCAL = 0, default).

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

  VU = V' * U;
  Ik = eye(k);

  beta = epsilon;
  delta = (1-epsilon)/2;
  B = Ik;
  D = -Ik/2;

  % A = epsilon*eye(n) + U*V';
  % X = A;
  % E = (eye(n) - A) / 2;

  % Xnew = beta*eye(n) + U*B*V';
  % Enew = delta*eye(n) + U*D*V';
  % norm(E - Enew) / norm(E)
  % norm(X - Xnew) / norm(X)

  if scal
    detA = abs(sqrt(det(epsilon*Ik + VU)))^(1/n);
    pre_mu = 1;
  end

  for iter = 1:maxit

    % disp('***');

    % Determinantal scaling.
    if scal
      mu = (sqrt(epsilon)/beta)^(1-k/n) *...
        detA / abs(det(beta*Ik + VU*B))^(1/n);
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

    delta_tilde = (delta + beta/2)/mu - (mu*beta)/2;
    D_tilde = (D + B/2)/mu - (mu*B)/2;

    beta_new = mu*beta + delta_tilde;
    B_new = mu*B + D_tilde;

    if norm(B-B_new, 1) / norm(B, 1) < tol
      break
    end

    S = D_tilde * VU;
    S(1:k+1:k^2) = S(1:k+1:k^2) + delta_tilde;

    T = VU * D_tilde;
    T(1:k+1:k^2) = T(1:k+1:k^2) + delta_tilde;

    % cond((beta_new*Ik + VU*B_new))
    % condB = cond(B_new)
    % condIns = cond(beta_new*Ik + VU*B_new)

    delta = -(delta_tilde^2) / (2*beta_new);
    % D = -((delta_tilde*Ik + S)*D_tilde -...
    %       S * B_new * ((beta_new*Ik + VU*B_new) \ T)) /...
    %     (2*beta_new);
    D = -((delta_tilde*Ik + S)*D_tilde -...
          S * B_new * (inv(beta_new*Ik + VU*B_new) * T)) /...
        (2*beta_new);

    beta = beta_new;
    B = B_new;

    % Etilde = (E+X/2)/mu - mu*X/2;
    % X = mu*X + Etilde;
    % E = -(Etilde * (X \ Etilde)) / 2;

    % Etildenew  = delta_tilde*eye(n) + U*D_tilde*V';
    % Xnew = beta*eye(n) + U*B*V';
    % Enew = delta*eye(n) + U*D*V';
    % Enew = -(1/2) * (delta_tilde*eye(n) + U*D_tilde*V') * ...
    %   ((beta*eye(n) + U*B*V') \ ...
    %   (delta_tilde*eye(n) + U*D_tilde*V'));

    % cond(beta*eye(n) + U*B*V')

    % Enew = -(1/2) * ((delta_tilde*eye(n) + U*D_tilde*V') / ...
    %                  (beta*eye(n) + U*B*V')) * ...
    %        (delta_tilde*eye(n) + U*D_tilde*V');
    % right = norm(E - Enew) / norm(E)

    % Enew = -(delta_tilde*eye(n) + U*D_tilde*V') * ...
    %        (eye(n) - U*B_new*((beta*Ik+VU)\V')) *...
    %        (delta_tilde*eye(n) + U*D_tilde*V') / (2*beta_new);
    % working_on = norm(E - Enew) / norm(E)

    % norm(Etilde - Etildenew) / norm(Etilde)
    % norm(X - Xnew) / norm(X)
    % norm(E - Enew) / norm(E)

  end

  X = U*B_new*V';
  X(1:n+1:n^2) = X(1:n+1:n^2) + beta_new;

end
