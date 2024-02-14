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
function T = profilex(frec, fmin, tau, testfeature)
    global FLAG_data;

    [np, ns, nr, maxfun] = size(frec);
    M = maxfun;

    T = NaN(np, ns, nr);
    f0 = -Inf(np, 1);

    for ip = 1:np
        f0(ip) = frec(ip, 1, 1, 1);
    end

    f_acc = zeros(np, ns);


    for ip = 1:np
        for is = 1:ns
            for ir = 1:nr
                if min(frec(ip, is, ir, 1:M)) <= tau * f0(ip) + (1 - tau) * fmin(ip) 
                    T(ip, is, ir) = find(frec(ip, is, ir, 1:M) <= tau * f0(ip) + (1 - tau) * fmin(ip), 1, 'first');
                    f_acc(ip, is) = (frec(ip, is, ir, T(ip, is, ir)) - f0(ip)) / (fmin(ip) - f0(ip));
                else
                    T(ip, is, ir) = NaN;
                end
            end
        end
    end

    sol = textread('solvers', '%s');
    for is = 1:ns
        sol{is} = regexprep(sol{is}, '_4test', '');
    end

    filename = 'problems1';
    fileID = fopen(filename);
    C = textscan(fileID, '%s %f');
    tstn = C{2};
    fclose(fileID);

    saveip = zeros(1000, 1);
    mm = 1;

    % mask = true(np, 1);
    % for ip = 1:np
    %     if (T(ip, 1) < T(ip, 5) - 2) || (T(ip, 2) <= T(ip, 5) - 2) || ...
    %        (T(ip, 3) < T(ip, 5) - 2) || (T(ip, 2) < T(ip, 4) - 15) || ...
    %        (T(ip, 1) < T(ip, 4) - 15)
    %         mask(ip) = false;
    %     else
    %         saveip(mm) = ip;
    %         mm = mm + 1;
    %     end
    % end
    % T = T(mask, :);

    save T;
    % save saveip;
    % save mm;
end
