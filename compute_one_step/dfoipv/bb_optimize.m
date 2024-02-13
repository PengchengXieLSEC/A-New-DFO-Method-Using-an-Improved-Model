function [res, iter, fhist] = bb_optimize(func, x_0, alg, clist, options)
  
  t1 = clock;
  if strcmp(lower(alg), 'dfo')
      [res, iter, fhist] = dfo_tr(func, x_0, clist, options);
  else
        res = minimize(func, x_0, alg, options);
        iter = [];
  end
  t2 = clock;
  time_consump = etime(t2, t1);
  % fprintf("Total time is %.3f seconds with %s method.",time_consump, alg);
end

