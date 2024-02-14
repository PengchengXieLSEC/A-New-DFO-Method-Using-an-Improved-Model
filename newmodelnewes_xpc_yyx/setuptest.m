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
function [xb, rhobeg, fopt, info] = setuptest (fun, n)
    % Set up x0, rhobeg, and fopt for the test.
    %character(len=*), intent(in) :: fun
    %integer(kind=4), intent(in) :: n
    %real(kind=8), intent(out) :: xb(n), rhobeg, fopt
    %integer(kind=4), intent(out) :: info

    %!      integer(kind=4), external :: mod
    %      character(len=len(fun)) :: ful
    %      integer(kind=4) :: i, j, k
    %      real(kind=8) :: ysrf(ceiling(sqrt((n))),+ceiling(sqrt((n)))), scaling, xo, pt, pi, theta
    xb = zeros(n, 1);
    xb(:) = 1.0e308;
    rhobeg = 1.0e0;
    fopt = 1.0e308;
    info = 0;

    ful = lower(fun);

    if strcmp(strtrim(ful), 'arglina')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'arglina4')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'arglinb')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'arglinc')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'argtrig')

        for i = 1:n
            xb(i) = 1.0e0 / (n);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'arwhead')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'bdqrtic')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'bdqrticp')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'bdvalue')

        for i = 1:n
            xb(i) = (i * (i - n - 1)) / ((n + 1) * (n + 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'brownal')
        xb(:)= 5.0e-1;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'broydn3d')
        xb(:)= -1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'broydn7d')
        xb(:)= -1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'brybnd')
        xb(:)= -1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'chainwoo')

        for i = 1:n

            if (i == 1 || i == 3)
                xb(i) = -3.0e0;
            elseif (i == 2 || i == 4)
                xb(i) = -1.0e0;
            else
                xb(i) = -2.0e0;
            end

        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'chebquad')

        for i = 1:n
            xb(i) = (i) / (n + 1);
        end

        rhobeg = 0.2e0 * xb(1);
    elseif strcmp(strtrim(ful), 'chpowellb')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'chpowells')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'chnrosnb')
        xb(:)= -1.0e0;
        % rhobeg = 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'chrosen')
        xb(:)= -1.0e0;
        rhobeg = 0.5e0;
    elseif strcmp(strtrim(ful), 'cosine')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'cragglvy')
        xb(:)= 2.0e0;
        xb(1) = 1.0e0;
        rhobeg = 1.0e0;
    elseif (strcmp(strtrim(ful), 'curly10') || strcmp(strtrim(ful), 'curly20') ...
            || strcmp(strtrim(ful), 'curly30'))

        for i = 1:n
            xb(i) = (1.0e-4) / ((n + 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'cube')
        xb(1) = -1.2e0;

        if (n >= 2)
            xb(2) = 1.0e0;
        end

        for i = 3:n
            xb(i) = xb(i - 2);
        end

        rhobeg = 1.0e0;
    elseif (strcmp(strtrim(ful), 'dixmaane') || strcmp(strtrim(ful), 'dixmaanf') ...
            || strcmp(strtrim(ful), 'dixmaang') || strcmp(strtrim(ful), 'dixmaanh') ...
            || strcmp(strtrim(ful), 'dixmaani') || strcmp(strtrim(ful), 'dixmaanj') ...
            || strcmp(strtrim(ful), 'dixmaank') || strcmp(strtrim(ful), 'dixmaanl') ...
            || strcmp(strtrim(ful), 'dixmaanm') || strcmp(strtrim(ful), 'dixmaann') ...
            || strcmp(strtrim(ful), 'dixmaano') || strcmp(strtrim(ful), 'dixmaanp'))

        for i = 1:n
            xb(i) = 2.0e0;
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'dqrtic')
        xb(:)= 2.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'edensch')
        xb(:)= 0.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'eg2')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'engval1')
        xb(:)= 2.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'errinros')
        xb(:)= -1.0e0;
        % rhobeg = 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'expsum')
        xb(:)= 0.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'extrosnb')
        xb(:)= -1.0e0;
        % rhobeg = 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'exttet')
        xb(:)= 1.0e-1;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'firose')
        xb(:)= -2.0e0;
        rhobeg = 0.5e0;
    elseif (strcmp(strtrim(ful), 'fletcbv2') || strcmp(strtrim(ful), 'fletcbv3'))

        for i = 1:n
            xb(i) = i / (n + 1);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'fletchcr')

        for i = 1:n
            xb(i) = 0.0e0;
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'fminsrf2')
        xb(:)= 0.0e0;
        ysrf = 0.0e0;
        k = floor(sqrt(n));

        if (k <= 1)
            info = 2;
        end

        for i = 1:k
            ysrf(1, i) = (4 * (i - 1)) / (k - 1) + 1.0e0;
            ysrf(i, 1) = (8 * (i - 1)) / (k - 1) + 1.0e0;
            ysrf(k, i) = (4 * (i - 1)) / (k - 1) + 1.0e0;
            ysrf(i, k) = (8 * (i - 1)) / (k - 1) + 1.0e0;
        end

        for i = 1:k

            for j = 1:k
                xb((i - 1) * k + j) = ysrf(i, j);
            end

        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'freuroth')
        xb(:)= 0.0e0;
        xb(1) = 0.5e0;

        if (n >= 2)
            xb(2) = -2.0e0;
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'genbrown')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'genhumps')
        xb(:)= -5.062e2;
        xb(1) = -5.060e2;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'genrose')

        for i = 1:n
            xb(i) = (i) / (n + 1);
        end

        % rhobeg = 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'indef')

        for i = 1:n
            xb(i) = (i) / (n + 1);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'integreq')

        for i = 1:n
            xb(i) = (i * (i - n - 1)) / ((n + 1) * (n + 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'liarwhd')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'lilifun3')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'lilifun4')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'morebv')

        for i = 1:n
            xb(i) = (i * (i - n - 1)) / ((n + 1) * (n + 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'morebvl')
        xb(:)= 1.0e0;
        % Suggested by Luksan, Matonoha, Vlcek, "Modified
        % CUTE probelems for unconstrained optimization.
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'ncb20')

        if (n < 20)
            xb(:)= 1.0e0;
        else
            xb(1:n - 10) = 1.0e0;
            xb(n - 9:n) = 1.0e0;
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'ncb20b')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'noncvxu2')

        for i = 1:n
            xb(i) = (i);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'noncvxun')

        for i = 1:n
            xb(i) = (i);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'nondia')
        xb(:)= -1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'nondquar')

        for i = 1:n

            if (mod(i, 2) == 1)
                xb(i) = 1.0e0;
            end

            if (mod(i, 2) == 0)
                xb(i) = -1.0e0;
            end

        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'penalty1')

        for i = 1:n
            xb(i) = (i);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'penalty2')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'penalty3')

        for i = 1:n
            xb(i) = (i) / (n + 1);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'penalty3p')
        xb(:)= 1.0e0;
        % Suggested by Powell, "The NEWUOA software
        % for unconstrained optimization without
        % derivatives", 2004, Section 8.
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'powellsg')

        for i = 1:n

            if (mod(i, 4) == 1)
                xb(i) = 3.0e0;
            end

            if (mod(i, 4) == 2)
                xb(i) = -1.0e0;
            end

            if (mod(i, 4) == 3)
                xb(i) = 0.0e0;
            end

            if (mod(i, 4) == 0)
                xb(i) = 1.0e0;
            end

        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'power')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'rosenbrock')
        xb(:)= -1.0e0;
        rhobeg = 0.5e0;
    elseif strcmp(strtrim(ful), 'sbrybnd')
        scaling = 1.2e1;

        for i = 1:n
            xb(i) = exp(scaling * (1 - i) / (n - 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'sbrybndl')
        scaling = 1.0e0;

        for i = 1:n
            xb(i) = exp(scaling * (1 - i) / (n - 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'schmvett')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'scosine')
        scaling = 1.0e0;

        for i = 1:n
            xb(i) = exp(scaling * (1 - i) / (n - 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'scosinel')
        scaling = 1.0e0;

        for i = 1:n
            xb(i) = exp(scaling * (1 - i) / (n - 1));
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'serose')
        xb(:)= -1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'sinquad')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'sparsine')
        xb(:)= 0.5e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'sparsqur')
        xb(:)= 0.5e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'spmsrtls')

        for i = 1:n
            xb(i) = 0.2e0 * sin((i)^2);
        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'sphrpts')

        if mod(n, 2) ~= 0
            info = 2;
        else

            for i = 1:n / 2
                xb(2 * i) = 0.0e0;
                xb(2 * i - 1) = 8.0e0 * acos(0.0e0) * (i) / (n);
            end

        end

        rhobeg = 1.0e0 / (n);
    elseif strcmp(strtrim(ful), 'srosenbr')

        for i = 1:n

            if (mod(i, 2) == 1)
                xb(i) = -1.2e0;
            end

            if (mod(i, 2) == 0)
                xb(i) = 1.0e0;
            end

        end

        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'stmod')
        xb(:)= 0.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'tointgss')
        xb(:)= 3.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'tointtrig')
        xb(:)= 1.0e0;
        rhobeg = 1.0e0;
    elseif strcmp(strtrim(ful), 'tquartic')
        xb(:)= 0.1e0;
        rhobeg = 1.0e0;
    elseif (strcmp(strtrim(ful), 'trigsabs') || strcmp(strtrim(ful), 'trigssqs'))

        for i = 1:n
            xo = pi * cos((n + i + 163) / 3.0e0);
            pt = pi * sin((n * i * 42)^(1.0e0/3.0e0));

            if strcmp(strtrim(ful), 'trigssqs')
                theta = (1.0e0 + sin((n * i * 42)^(1.0e0/3.0e0))) / 1.0e0;
                theta = 1.0e1^(-theta);
                xb(i) = (xo + pt * 1.0e-1) / theta;
            else
                xb(i) = xo + pt * 1.0e-1;
            end

        end

        rhobeg = 1.0e-1;
    elseif strcmp(strtrim(ful), 'trirose1')
        xb(:)= -1.0e0;
        rhobeg = 0.5e0;
    elseif strcmp(strtrim(ful), 'trirose2')
        xb(:)= 12.0e0;
        rhobeg = 2.0e0;
    elseif strcmp(strtrim(ful), 'vardim')

        for i = 1:n
            xb(i) = (n - i) / (n);
        end

        rhobeg = 1.0e0 / (2 * n);
    elseif strcmp(strtrim(ful), 'woods')

        for i = 1:n

            if (mod(i, 2) == 1)
                xb(i) = -3.0e0;
            end

            if (mod(i, 2) == 0)
                xb(i) = -1.0e0;
            end

        end

        rhobeg = 1.0e0;
    else
        info = 1;
    end

end
