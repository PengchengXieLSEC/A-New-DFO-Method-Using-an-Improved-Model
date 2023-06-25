%Codes for the paper entitled "A New Derivative-free Method Using an Improved Under-determined Quadratic Interpolation Model"
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn


function [coeff, s, nrms, w] = trust_sub_compute_step(alpha, eigval, coeff, V, lam)
    w = eigval + lam;
    arg1 = (w == 0 & alpha == 0);
    arg2 = (w == 0 & alpha ~= 0);
    coeff(w ~= 0) = alpha(w ~= 0) ./ w(w ~= 0);
    coeff(arg1) = 0;
    coeff(arg2) = inf;
    coeff(isnan(coeff)) = 0;
    s = V * coeff;
    nrms = norm(s);
end

