function X = sqrtm_direct(alpha, U, V)
% Compute the square root of ALPHA*I + U*V' using formula (1.5) in our
  [n,k] = size(U);
  M = V'*U;
  X = M;
  X(1:k+1:k^2) = X(1:k+1:k^2) + alpha;
  X = sqrtm(X);
  X(1:k+1:k^2) = X(1:k+1:k^2) + sqrt(alpha);
  X = U * (X \ V');
  X(1:n+1:n^2) = X(1:n+1:n^2) + sqrt(alpha);
end
