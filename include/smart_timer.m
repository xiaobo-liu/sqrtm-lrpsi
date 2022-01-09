function [X, time] = smart_timer(f, varargin)
% SMART_TIMER Measure execution time of function.
%    [X,TIME]=SMART_TIMER(F,VARARGIN) returns the average execution time of the
%    execution of the function F evaluated at VARARGIN. The function is executed
%    only once if the first run takes more than second, and multiple times
%    otherwise.
  init = tic;
  X = f(varargin{:});
  time = toc(init);
  if time < 1
    timefun = @()(f(varargin{:}));
    time = timeit(timefun);
  end
end