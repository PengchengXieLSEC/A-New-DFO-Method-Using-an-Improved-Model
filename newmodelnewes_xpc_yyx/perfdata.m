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
function perfdata(feature)

    options = [];
    test_options = [];

    frec = testalg('ALL', 'ALL', test_options);
    sol = textread('solvers', '%s');
    ns = length(sol);
    global prob
    np = length(prob);
    History_Fs = cell(1, ns);
    
    for is = 1:ns
        History_Fs{is} = readmatrix(['history', num2str(is), '.csv']);
    end
    
    NFs = readmatrix('NF_data.csv');
    
    for is = 1:ns
        FLAG = 1;
        NF = NFs(:, is);
        History_F = History_Fs{is};
        
        for ip = 1:np
            frec(ip, is, 1, 1:NF(ip)) = History_F(FLAG:FLAG + NF(ip) - 1);
            FLAG = FLAG + NF(ip);
            
            if NF(ip) == 110000
                continue
            else
                frec(ip, is, 1, NF(ip):110000) = History_F(FLAG - 1);
            end
        end
    end

    save frec

    [np, ns, ~, ~] = size(frec);
    fmin = NaN(np, 1);

    for ip = 1:np
        fmin(ip) = min(min(min(frec(ip, :, 1, :))));
    end


    for tau = 10 .^ (-5:-1:-5)
        profilex(frec, fmin, tau, 'plain');
    end

    return;
end
