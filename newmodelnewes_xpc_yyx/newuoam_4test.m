function [x, f, exitflag, output] = newuoam_4test(fun, x0, options, n)

  if (~exist('newuoam', 'dir'))
    warning('newuoam seems unavailable on your computer. Forget about testing it.');
    x = x0;
    f = fun(x0);
  else
    cd 'newuoam'
    [x, nf, f] = newuoam(fun, n, 2 * n + 1, x0, options.rhobeg, options.rhoend, 0, options.maxfun);
    cd '..'
    exitflag = 0;
    output.funcCount = nf;
  end
end
