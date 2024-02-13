function [lambda_0] = quad_Frob_compute_coeffs(W, tol_svd, b, option)
%QUAD_FROB_COMPUTE_COEFFS 此处显示有关此函数的摘要
%   此处显示详细说明
    if strcmp(option, 'partial')
       [U, S, VT] = svd(W);
    end
    
    if cond(W) < 1.0e+10
        lambda_0 = W \ b;
    else
%     Make sure the condition number is not too high
       indices = (S < tol_svd);
       S(indices) = tol_svd;
       Sinv = zeros(size(S,1), size(S,1));
       Sinv(~indices) = 1./S(~indices);
%         Get the coefficients
       lambda_0 = VT * Sinv * U' * b;
    end
end

