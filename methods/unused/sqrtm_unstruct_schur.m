function X = sqrtm_unstruct_schur(epsilon, U, V)
% Compute the square root of EPSILON*I + U*V' directly using the Schur 
% method.
  n = size(U, 1);
  A = epsilon*eye(n) + U*V';
  X = sqrtm(A);
end
