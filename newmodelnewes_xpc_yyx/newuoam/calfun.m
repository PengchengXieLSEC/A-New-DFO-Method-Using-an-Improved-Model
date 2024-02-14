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
function [F] = calfun (N, X)
  Y = zeros(10);
  for J = 1:N
    Y(1, J) = 1.0e0;
    Y(2, J) = 2.0e0 * X(J) - 1.0e0;
  end
  for I = 2:N
    for J = 1:N
      Y(I + 1, J) = 2.0e0 * Y(2, J) * Y(I, J) - Y(I - 1, J);
    end
  end
  F = 0.0e0;
  NP = N + 1;
  IW = 1;
  for I = 1:NP
    SUM = 0.0e0;
    for J = 1:N
      SUM = SUM + Y(I, J);
    end
    SUM = SUM / (N);
    if (IW > 0)
      SUM = SUM + 1.0e0 / (I * I - 2 * I);
    end
    IW = -IW;
    F = F + SUM * SUM;
  end
end
