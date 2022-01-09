function X = sqrtm_unstruct(epsilon, U, V)
% Compute the square root of EPSILON*I + U*V' using diagonalization.
  n = size(U, 1);
  A = epsilon*eye(n) + U*V';
  if isequal(U, V) % symmetric case
    X = sqrtm(A);
    % [P, d] = eig(A, 'vector');
    % X = P * (sqrt(max(d,0)) .* P');
    % X = (X + X')/2;
  else
    X = sqrtm(A);
  end
end
