function X = sqrtm_dir_db_prod(alpha, U, V)
% Compute the square root of ALPHA*I + U*V' using formula in Cor. 1.3,
% where the k-by-k matrix square root is computed by the product form of
% the DP iteration.
  [n,k] = size(U);
  M = V'*U;
  X = M;
  X(1:k+1:k^2) = X(1:k+1:k^2) + alpha;
  X = sqrtm_dbp(X);
  X(1:k+1:k^2) = X(1:k+1:k^2) + sqrt(alpha);
  X = U * (X \ V');
  X(1:n+1:n^2) = X(1:n+1:n^2) + sqrt(alpha);
end
