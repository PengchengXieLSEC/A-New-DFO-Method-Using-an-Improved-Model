function [res] = get_results(func, x_0, alg, options, para)
 
    [res] = bb_optimize(func, x_0, alg, options, para);

    fprintf ("\n"+"Printing result for function " + func + ": \n")
    fprintf ("best point: %.4f, with obj: %.2f  ", res.x', res.fun)
    fprintf ("\n-------------" + alg + " Finished ----------------------\n")
end

