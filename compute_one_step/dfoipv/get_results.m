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
function [res] = get_results(func, x_0, alg, options)
    c6 = [[1/3, 1/3, 1/3];...
        [1, 0, 0]; ...
        [0, 1, 0]; ...
        [0, 0, 1]; ...
        [1/2, 1/2, 0]; ...
        [0, 1/2, 1/2]];
    length = size(c6,1);
    for j = 1:length
        [res, it, fhist] = bb_optimize(func, x_0, "DFO", c6(j, 1:end), options);
        %[res] = bb_optimize(func, x_0, alg, options);

    % fprintf ("\n"+"Printing result for function " + func + ": \n")
    % fprintf ("best point: %.4f, with obj: %.2f  ", res.x', res.fun)
    % fprintf ("\n-------------" + alg + " Finished ----------------------\n")
end

