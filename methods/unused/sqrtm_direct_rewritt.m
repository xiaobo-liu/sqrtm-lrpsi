function X = sqrtm_direct_rewritt(alpha, U, V)
% Compute the square root of EPSILON*I + U*V' using the direct formula in
% Theorem 1.35. The formula is rewritten to mitigate severe cancellation
% for small U*V'.
  [n,k] = size(U);
  M = V'*U;
  X = M;
  X(1:k+1:k^2) = X(1:k+1:k^2) + alpha;
  X = sqrtm(X);
  X(1:k+1:k^2) = X(1:k+1:k^2) + sqrt(alpha);
  X = U * ((X * M) \ M) * V';
  X(1:n+1:n^2) = X(1:n+1:n^2) + sqrt(alpha);
end