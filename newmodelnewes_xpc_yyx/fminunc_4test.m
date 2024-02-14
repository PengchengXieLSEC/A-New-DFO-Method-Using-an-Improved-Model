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
function [x, f, exitflag, output] = fminunc_4test(func, x0, options, ~)

  tol = options.rhoend;
 % maxfun = options.maxfun;
h = sqrt(eps);

fminunc_options = optimoptions(@fminunc, 'FiniteDifferenceStepSize', 0.1, 'Algorithm', 'quasi-newton', 'FunctionTolerance', tol, 'StepTolerance', tol, 'OptimalityTolerance', tol, 'MaxFunctionEvaluations', 100*length(x0), 'MaxIterations', 100*length(x0), 'Display', 'none');

[x, f, exitflag, output] = fminunc(func, x0, fminunc_options);

end
