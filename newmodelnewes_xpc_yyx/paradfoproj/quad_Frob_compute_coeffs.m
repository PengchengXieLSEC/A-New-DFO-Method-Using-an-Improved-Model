% A DFO Method Using a New Model
% Codes for the paper entitled
% "A Derivative-free Method Using a New Under-determined Quadratic Interpolation Model"
% Copyright: Pengcheng Xie & Ya-xiang Yuan

% Connect: xpc@lsec.cc.ac.cn
% A DFO Method Using a New Model
% ----------------------------------------------------------
% License Information

% ----------------------------------------------------------
% This code is distributed under the MIT License.
% You should have received a copy of the MIT License along
% with this program. If not, see <https://opensource.org/licenses/MIT>.

% For further information or questions, contact the authors
% via the provided email address.
% ----------------------------------------------------------
% Code Version Information

% ----------------------------------------------------------
% Version: 1.0
% Changes: Initial release.
% ----------------------------------------------------------

% ----------------------------------------------------------
% References
% ----------------------------------------------------------
% For more information, refer to the paper:

% "A Derivative-free Method Using a New Under-determined Quadratic Interpolation Model"
% by Pengcheng Xie & Ya-xiang Yuan.
%
% If you use this code in your research, please cite the above paper.

% ----------------------------------------------------------
% ----------------------------------------------------------
% Contributors
% ----------------------------------------------------------

% This code was written by Pengcheng Xie & Ya-xiang Yuan.
% ----------------------------------------------------------
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

