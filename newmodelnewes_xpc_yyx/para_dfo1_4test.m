function [x, f, alg, func] = para_dfo1_4test(func, x_0, ~, ~)
  if (~exist('paradfo', 'dir'))
    warning('paradfo seems unavailable on your computer. Forget about testing it.');
    x = x_0;
    f = fun(x_0);
  else
    options = struct("maxfev", 100*length(x_0), "init_delta", 1, "tol_delta", 1e-6, "tol_f", 1e-6, "tol_norm_g", 1e-6, "sample_gen", "auto", "verbosity", 0);
    alg = "DFO";
    cd 'paradfo'
    [res, ~, f] = bb_optimize(func, x_0, alg, options);
    cd '..'
    x = res.x;
  end
end
