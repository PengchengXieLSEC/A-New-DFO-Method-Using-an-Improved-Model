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
function [x, f, exitflag, output] = newuoam_4test(fun, x0, options, n)

  if (~exist('newuoam', 'dir'))
    warning('newuoam seems unavailable on your computer. Forget about testing it.');
    x = x0;
    f = fun(x0);
  else
    cd 'newuoam'
    [x, nf, f] = newuoam(fun, n, 2 * n + 1, x0, options.rhobeg, options.rhoend, 0, options.maxfun);
    cd '..'
    exitflag = 0;
    output.funcCount = nf;
  end
end
