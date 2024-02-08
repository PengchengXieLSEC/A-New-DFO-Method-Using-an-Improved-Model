function [x, f, exitflag, output] = fminsearch_4test(func, x0, options, ~)

  tol = options.rhoend;
  fminsearch_options = optimset('MaxFunEvals', 100*length(x0), 'TolFun', tol, 'TolX', tol, 'Display', 'none');
  [x, f, exitflag, output] = fminsearch(func, x0, fminsearch_options);

end
