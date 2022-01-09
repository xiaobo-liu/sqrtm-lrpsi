function s = latex_sci_not(x, fpdigits, expdigits)
% LATEX_SCI_NOT Format real number as a LaTeX string in scientific notation.
%    S=LATEX_SCI_NOT(X,FPDIGITS,EXPDIGITS) returns a string that represents the
%    integer X in scientific notation. The number has FPDIGITS decimal digits
%    after the decimal point, and the exponenet is repersented using EXPDIGITS
%    decimal digits. If the exponent cannot be represented using EXPDIGITS
%    decimal digits a warning is emitted. The default value of the arguments
%    FPDIGITS and EXPDIGITS are 1 and 2, respectively.

  if nargin < 2
    fpdigits = 1;
  end
  if nargin < 3
    expdigits = 2;
  end

  exp10 = floor(log10(abs(x)));
  xround = round(x*(10^(-exp10+fpdigits))) * 10^(exp10-fpdigits);
  exp10 = floor(log10(abs(xround)));
  if exp10 == 0
    currexpdigits = 1;
  else
    currexpdigits = floor(log10(abs(exp10)))+1;
  end

  if currexpdigits < expdigits % Add phantom trailing zeros.
    s = sprintf('$%%.%df\\\\times 10^{%%d\\\\phantom{%%0%dd}}$',...
                fpdigits, expdigits-currexpdigits);
    s = sprintf(s, xround*(10^-exp10), exp10, 0);
  else
    if currexpdigits > expdigits
      warning('More than %d digits are required to typeset the exponent of %d.',...
              expdigits, x);
    end
    s = sprintf('$%%.%df\\\\times 10^{%%d}$', fpdigits);
    s = sprintf(s, xround*(10^-exp10), exp10);
  end
end