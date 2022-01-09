function err = compute_rel_err(Xref, X, mynorm, use_high_prec)
% COMPUTE_REL_ERR   Compute relative error in the chosen matrix norm.

%    ERR=COMPUTE_REL_ERR(XREF,X,MYNORM,USE_HIGH_PREC) computes the relative
%    error
%        || X - XREF ||_MYNORM / || XREF ||_MYNORM.
%    If USE_HIGH_PREC is set to TRUE the matrix subtraction is performed in
%    higher precision, otherwise it is computed in the same precision as X. The
%    parameters MYNORM and USE_HIGH_PREC default to 2 and TRUE, respectively.

  if nargin < 2 || nargin > 3
    error('This function requires two or three arguments.');
  end
  if nargin < 3
    mynorm = 2;
  end
  if nargin < 4
    use_high_prec = true;
  end

  if use_high_prec
    ndigits = mp.Digits;
    if isa(X, 'single')
      highprec = @(x)(double(x));
    elseif isa(X, 'double')
      mp.Digits(34);
      highprec = @(x)(mp(x));
    elseif isa(X, 'double')
      mp.Digits(2*ndigits);
      highprec = @(x)(mp(x));
    end
    err = norm(highprec(X) - highprec(Xref), mynorm) /...
          norm(highprec(Xref), mynorm);
    mp.Digits(ndigits);
  else
    err = norm(X - Xref, mynorm) / norm(Xref);
  end
end
