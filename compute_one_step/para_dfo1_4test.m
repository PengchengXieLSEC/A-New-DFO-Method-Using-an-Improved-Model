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
function [x, f, alg, func] = para_dfo1_4test(func, x_0, ~, ~)
  if (~exist('paradfo', 'dir'))
    warning('paradfo seems unavailable on your computer. Forget about testing it.');
    x = x_0;
    f = fun(x_0);
  else
    options = struct("maxfev", 100*length(x_0), "init_delta", 1, "tol_delta", 1e-6, "tol_f", 1e-6, "tol_norm_g", 1e-6, "sample_gen", "manual", "verbosity", 0);
    alg = "DFO";
    cd 'paradfo'
    [res, ~, f] = bb_optimize(func, x_0, alg, options);
    cd '..'
    x = res.x;
  end
end
