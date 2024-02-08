% separate the functions based on their dimension. This is merely
% done to ensure the starting point x_0 will later have the
% correct dimension
%nd_func_names = [sphere, rosen]  # functions from R^n -> R
clc
clear all

nd_func_names = {'ackley'};
%td_func_names = [ackley, booth, bukin, beale]   # functions from R^2 -> R
td_func_names = {'ackley'};

all_func_names = [td_func_names, nd_func_names];

% Run all the algorithms and problems with given starting points
% Specify the starting point and options. For example, try the following
% options.
for i = size(all_func_names,2)
    func = string(all_func_names(i));
    if ismember(func, td_func_names)
        x_0 = [1.3, 0.7];
    else
        x_0 = [1.3, 0.7, 0.8, 1.9, 1.2];
    end
    % fprintf("\n\n********* Function " + func + "********")
    alg = "DFO";
    options = struct("maxfev", 100, "init_delta", 1, "tol_delta", 1e-6, "tol_f", 1e-6, "tol_norm_g", 1e-5, "sample_gen", "auto", "verbosity", 1);
    % 访问特定元素用options.maxfev
    res = get_results(func, x_0, alg, options);
end