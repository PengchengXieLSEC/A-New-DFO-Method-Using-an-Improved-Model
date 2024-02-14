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
function frec = testalg(solver, problem, options)

    global fvals fval_history
    global F_data X_data NF_data prob ip is



    if strcmpi('ALL', problem)
        filename = 'problems1';
        fileID = fopen(filename);
        C = textscan(fileID, '%s %f');
        tstn = C{2};
        prob = C{1};
        fclose(fileID);
    else
        prob = {problem};
    end

    rhoend = 1e-6;
    tol = rhoend;

    nr = 1;


    np = length(prob);

    if strcmpi('ALL', solver)
        sol = textread('solvers', '%s');
    else
        sol = {solver};
    end
    ns = length(sol);

    F_data = zeros(np, ns);
    X_data = cell(np, 1);
    for iip = 1:np
        fvals{iip} = zeros(1, ns * tstn(ip));
    end
    NF_data = zeros(np, ns);

    frec = NaN(np, ns, nr, 60000);
    fvals = cell(1, ns);
    for iis = 1:ns
        fvals{iis} = [];
    end

    for ip = 1:np
        disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp(strcat(int2str(ip), '. ', prob{ip}, ':'));

            maxfun = 100 * tstn(ip); 

        [x0, rhobeg, ~, ~] = setuptest(prob{ip}, tstn(ip));
        for is = 1:ns
            solv = sol{is};
            if np * nr <= 20
                disp('-------------------------------------------------------');
            end
            display(strcat(sol{is}, ':'));

            for ir = 1:nr
                if exist('solv_options', 'var') == 1
                    clear solv_options;
                end
                solv_options.rhobeg = rhobeg;
                solv_options.rhoend = rhoend;
                solv_options.tol = tol;
                solv_options.maxfun = maxfun;

                fval_history = [];
                if strcmp(solv, 'fminunc_4test') || strcmp(solv, 'cmaes_4test') || ...
                   strcmp(solv, 'newuoam_4test') || strcmp(solv, 'fminsearch_4test') || ...
                   strcmp(solv, 'nms_4test') || strcmp(solv, 'dfom_4test') || ...
                   strcmp(solv, 'dfoc_4test') || strcmp(solv, 'dfo1_4test') || ...
                   strcmp(solv, 'para_dfo1_4test') || strcmp(solv, 'para_dfo2_4test')

                    Func_wraper = @(x)(evalobjfun(prob{ip}, x, tstn(ip)));
                    solverHub = str2func(solv);
                    [x, ~, ~, ~] = solverHub(Func_wraper, x0, solv_options, tstn(ip));
                end

                nf = length(fval_history);

                X_data{ip}(is * tstn(ip) - (tstn(ip) - 1):is * tstn(ip)) = x';
                F_data(ip, is) = evalfun(prob{ip}, x, tstn(ip));
                NF_data(ip, is) = nf;

                nf = min(nf, maxfun);

                frec(ip, is, ir, 1:nf) = fval_history(1:nf);
                frec(ip, is, ir, nf + 1:maxfun) = fval_history(nf);
                if np * nr <= 20
                    fprintf('nf = %d\n', nf);
                    minfval = min(fval_history);
                    fprintf('minfval = %.6E\n', minfval);
                    if tstn(ip) <= 30
                        xx = sprintf(' %.3E ', x');
                        fprintf('xbest = %s\n', xx);
                    end
                end
            end
        end
        disp('--------------------------------------------------------------------');
        disp('--------------------------------------------------------------------');
        disp(['Solving problem' num2str(ip) '  Process：' num2str(ip / np * 100) '%']);
        disp('--------------------------------------------------------------------');
        disp('--------------------------------------------------------------------');
    end

    writematrix(F_data, 'F_data.csv');
    writecell(X_data, 'X_data.csv');
    writematrix(NF_data, 'NF_data.csv');
    for iis = 1:ns
        writematrix(fvals{iis}, ['history', num2str(iis), '.csv']);
    end

end

function f = evalobjfun(fun, x, n)
    global fvals is fval_history;

    [f, ~] = evalfun(fun, x, n);
    fval_history = [fval_history, f];
    fvals{is} = [fvals{is}; f];
end
