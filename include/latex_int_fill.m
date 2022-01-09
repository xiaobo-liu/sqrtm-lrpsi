function s = latex_int_fill(x, ndigits)
% LATEX_INT_FILL   Format integer as a LaTeX string.
%    S=LATEX_INT_FILL(X,NDIGITS) returns a string that represents the integer X
%    using NDIGITS digits in math mode LaTeX. If X has more than NDIGITS digits,
%    the string representing X is returned and a warning is emitted.

  if nargin < 2 || nargin > 2
    error('Exctly two arguments must be specified.');
  end
  if ~(isscalar(x) && isreal(x) && round(x) == x)
    error('The first argument must be an integer.')
  end

  ndigitsx = floor(log10(x))+1;
  if ndigitsx < ndigits % Fill with phanotm zeros.
    s = sprintf('$\\\\phantom{%%0%dd}%%d$',...
                ndigits - ndigitsx);
    s = sprintf(s, 0, x);
  else
    if ndigitsx > ndigits
      warning('More than %d digits are required to typeset %d.',...
              ndigits, x);
    end
    s = sprintf('$%d$', x);
  end
end