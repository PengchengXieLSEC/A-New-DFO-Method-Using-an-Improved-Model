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
function [output] = funcs_def(func_name,x)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
  if strcmp(func_name, 'rosen')
      output = sum(100.0 * (x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
  elseif  strcmp(func_name, 'sphere')
      output = sum(x.^2);
  elseif  strcmp(func_name, 'ackley')
      e = exp(1);
      output = -20 * exp(-0.2 * sqrt(0.5 * sum(x.^2))) - ( - exp(0.5 * (sum(cos(2*pi*x))))) + e + 20;
  elseif  strcmp(func_name, 'beale')
      x1 = x(1);
      x2 = x(2);
      output = (1.5 - x1 + x1*x2).^2 + (2.25 - x1 + x1*x2.^2).^2 + ((2.625 - x1 + x1*x2.^3).^2);
  elseif  strcmp(func_name, 'booth')
      x1 = x(1);
      x2 = x(2);
      output = ((x1+2*x2-7).^2 + (2*x1+x2-5).^2);
  elseif  strcmp(func_name, 'bukin')
      x1 = x(1);
      x2 = x(2);
      output = (100 * sqrt(abs(x2 - 0.01*x1.^2)) + 0.01 * (abs(x1+10)));
  elseif  strcmp(func_name, 'norm2')
      x1 = x(1);
      x2 = x(2);
      output = x1^2 + x2^2;
  end
end

