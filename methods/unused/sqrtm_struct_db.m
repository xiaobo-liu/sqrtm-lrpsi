function [Y, iter] = sqrtm_struct_db(epsilon, U, V, scal, tol, maxit)
% Compute the square root of EPSILON*I + U*V' using a structured version of the
% Denman-Beavers iteration. SCAL specifies whether scaling is used (SCAL = 1) or
% not (SCAL = 0, default).

% TODO: Can the cost of the first few iterations be reduced in this case?

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
  gamma = epsilon;
  C = eye(k);
  delta = 1;
  D = zeros(k);

  % A = epsilon*eye(n) + U*V';
  % Y = gamma*eye(n) + U*C*V';
  % Z = delta*eye(n) + U*D*V';

  % norm(Y - A)
  % norm(Z - eye(n))

  if scal
    detA = abs(sqrt( det(epsilon*eye(k) + M) ))^(1/n);
    pre_mu = 1;
  end


  for iter = 1:maxit

    MC = M*C;

    % Determinantal scaling.
    if scal
      mu = (sqrt(epsilon)/gamma)^(1-k/n) * detA / abs(det(gamma*eye(k) + MC))^(1/n);
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

    % if mu ~= 1
    %  disp('!')
    % end
    % Ynew = (mu*Y+inv(Z)/mu)/2;
    % Z = (mu*Z+inv(Y)/mu)/2;
    % Y = Ynew;

    mugamma = mu*gamma;
    mudelta = mu*delta;

    gamma_new = (mugamma + 1/mudelta)/2;
    C_new = (mu/2)*C - (D/(delta*eye(k) + M*D))/(2*mudelta);

    delta_new = (mudelta + 1/mugamma)/2;
    D_new = (mu/2)*D - (C/(gamma*eye(k) + MC))/(2*mugamma);

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
    % X = gamma*eye(n) + U*C*V';
    % norm(Y^2-A) / norm(A)

  end

  % Y  = gamma*eye(n) + U*C*V';
  Y = U*C_new*V';
  Y(1:n+1:n^2) = Y(1:n+1:n^2) + gamma_new;

end
