%Codes for the paper entitled "A New Derivative-free Method Using an Improved Under-determined Quadratic Interpolation Model"
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function [res, iter, fhist] = bb_optimize(func, x_0, alg, options, para)
  
  t1 = clock;
  if strcmp(lower(alg), 'dfo')
      [res, iter, fhist] = dfo_tr(func, x_0, options, para);
  else
        res = minimize(func, x_0, alg, options);
        iter = [];
  end
  t2 = clock;
  time_consump = etime(t2, t1);
  fprintf("Total time is %.3f seconds with %s method.",time_consump, alg);
end

