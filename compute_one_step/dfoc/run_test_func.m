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
    fprintf("\n\n********* Function " + func + "********")
    alg = "DFO";
    options = struct("maxfev", 100, "init_delta", 1, "tol_delta", 1e-6, "tol_f", 1e-6, "tol_norm_g", 1e-5, "sample_gen", "auto", "verbosity", 1);
    % 访问特定元素用options.maxfev
    out = zeros(4096,12);% 1~6为para，7、8为x4，9、10为x5，11、12为G+omega是否正定
    % 第二小的数不等于0
%     para = [c_1, c_2, omega_1, omega_2 , mu, theta];
    para = [1, 1, 1, 1, 1, 1];
    [res] = get_results(func, x_0, alg, options, para);
end