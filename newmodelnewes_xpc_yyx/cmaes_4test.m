function [x, f, alg, func] = cmaes_4test(func, x_0, ~, ~)
  if (~exist('CMA-ES', 'dir'))
    warning('CMA-ES seems unavailable on your computer. Forget about testing it.');
    x = x_0;
    f = func(x_0);
  else
    alg = "CMA-ES";
    cd 'CMA-ES'
    [x, f] = cmaes(func, x_0, 1);
    cd '..'
  end
end
