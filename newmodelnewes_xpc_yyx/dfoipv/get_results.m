function [res] = get_results(func, x_0, alg, options)
    c6 = [[1/3, 1/3, 1/3];...
        [1, 0, 0]; ...
        [0, 1, 0]; ...
        [0, 0, 1]; ...
        [1/2, 1/2, 0]; ...
        [0, 1/2, 1/2]];
    length = size(c6,1);
    for j = 1:length
        [res, it, fhist] = bb_optimize(func, x_0, "DFO", c6(j, 1:end), options);
        %[res] = bb_optimize(func, x_0, alg, options);

    % fprintf ("\n"+"Printing result for function " + func + ": \n")
    % fprintf ("best point: %.4f, with obj: %.2f  ", res.x', res.fun)
    % fprintf ("\n-------------" + alg + " Finished ----------------------\n")
end

