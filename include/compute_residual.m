function res = compute_residual(A, X, mynorm, use_high_prec)
% COMPUTE_RESIDUAL   Compute residual of matrix square root.
%    RES=COMPUTE_RESIDUAL(A,X,MYNORM,USE_HIGH_PREC) computes the residual
%        || X^2 - A ||_MYNORM / || A ||_MYNORM.
%    If CHANGE_PREC is set to TRUE the matrix product is computed in higher
%    precision, otherwise it is computed in the same precision as X. The
%    parameters MYNORM and USE_HIGH_PREC default to 2 and TRUE,
%    respectively.

  if nargin < 2 || nargin > 4
    error('This function requires two or three arguments.');
  end
  if nargin < 3
    mynorm = 2; % now use 2-norm by default
  end
  if nargin < 4
    use_high_prec = true;
  end

  % Higher precision used for computing X^2.
  if use_high_prec
    ndigits = mp.Digits;
    if isa(X, 'single')
      highprec = @(x)(double(x));
    elseif isa(X, 'double')
      mp.Digits(34);
      highprec = @(x)(mp(x));
    elseif isa(X, 'mp')
      mp.Digits(2*ndigits);
      highprec = @(x)(mp(x));
    end
    % highprec = @(x)(x);
    res = cast(norm(highprec(X)*highprec(X) - highprec(A), mynorm) /...
               norm(highprec(A), mynorm), class(X));
    mp.Digits(ndigits);
  else
    res = norm(X*X - A, mynorm) / norm(A, mynorm);
  end
end
