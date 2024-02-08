function [x, f, exitflag, output] = fminunc_4test(func, x0, options, ~)

  tol = options.rhoend;
 % maxfun = options.maxfun;
h = sqrt(eps);

fminunc_options = optimoptions(@fminunc, 'FiniteDifferenceStepSize', 0.1, 'Algorithm', 'quasi-newton', 'FunctionTolerance', tol, 'StepTolerance', tol, 'OptimalityTolerance', tol, 'MaxFunctionEvaluations', 100*length(x0), 'MaxIterations', 100*length(x0), 'Display', 'none');

[x, f, exitflag, output] = fminunc(func, x0, fminunc_options);

end
