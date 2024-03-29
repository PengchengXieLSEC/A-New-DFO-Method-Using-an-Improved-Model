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
function [value] = trust_sub_secular_eqn(lambda_0, eigval, alpha, delta)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    m = size(lambda_0,1);
    n = size(eigval,1);
    unn = ones(n, 1);
    unm = ones(m, 1);
    M = eigval * unm' + unn * lambda_0';
    MC = M;
    MM = alpha * unm';
    M(M ~= 0.0) = MM(M ~= 0.0) ./ M(M ~= 0.0);
    M(MC == 0.0) = inf;
    M = M.*M;
    value = sqrt(unm / (M' * unn));

    if size(value(value == inf),1)
        inf_arg = (value == inf);
        value(inf_arg) = zeros(size(value(inf_arg),1), 1);
    end

    value = (1.0./delta) * unm - value;
end

