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
function [fun, info] = evalfun(name, x, n)
    %% TestProblemF

    INFINITY = 1.0e308;
    f = INFINITY;
    info = 0; % info = 0 means successful evaluation
    y = zeros(n + 2, 1);
    y(0 + 1) = 0.0e0;
    y(n + 1 + 1) = 0.0e0;

    for i = 1:n
        y(i + 1) = x(i);
    end

    ychebq = zeros(n + 1, n);
    ysrf = zeros(ceil(sqrt(n)), ceil(sqrt(n)));
    ytmp = zeros(n, 1);
    par = zeros(n, 1);
    t = zeros(n, 1);
    MST = zeros(n, n);
    MaS = zeros(2 * n, n);
    MaC = zeros(2 * n, n);
    sinx = zeros(n, 1);
    cosx = zeros(n, 1);
    so = zeros(n, 1);
    co = zeros(n, 1);
    d = zeros(2 * n, 1);
    f1 = 0;
    f2 = 0;
    f3 = 0;
    f4 = 0;
    f5 = 0;
    name = lower(name);
    %%
    if (strcmp(name, 'arglina'))
        %% see j. more, 1981, testing unconstrained optimization software,;
        % problem 32, also cuter.;
        m = fix(2 .* n);
        tmp = 0.0e0;

        for i = 1:n
            tmp = tmp + y(i + 1);
        end

        f = 0.0e0;

        for i = 1:n
            f = f + (y(i + 1) - 2.0e0 ./ real(m) .* tmp - 1.0e0).^2;
        end

        for i = n + 1:m
            f = f + (-2.0e0 ./ real(m) .* tmp - 1.0e0).^2;
        end

    elseif (strcmp(name, 'arglina4'))
        %!!!!! DIFFERENT FROM ARGLINA. Quartic instead of  Quadratic. !!!!
        %%
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 32, also cuter.;
        m = fix(2 .* n);
        tmp = 0.0e0;

        for i = 1:n
            tmp = tmp + y(i + 1);
        end

        f = 0.0e0;

        for i = 1:n
            f = f + (y(i + 1) - 2.0e0 ./ real(m) .* tmp - 1.0e0).^4;
        end

        for i = n + 1:m
            f = f + (-2.0e0 ./ real(m) .* tmp - 1.0e0).^4;
        end

    elseif (strcmp(name, 'arglinb'))
        %% arglinb
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 33, also cuter.;
        m = fix(2 .* n);
        tmp = 0.0e0;

        for i = 1:n
            tmp = tmp + real(i) .* y(i + 1);
        end

        f = 0.0e0;

        for i = 1:m
            f = f + (real(i) .* tmp - 1.0e0).^2;
        end

    elseif (strcmp(name, 'arglinc'))
        %% arglinc
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 34, also cuter.;
        m = fix(2 .* n);
        tmp = 0.0e0;

        for i = 2:n - 1
            tmp = tmp + real(i) .* y(i + 1);
        end

        f = 1.0e0;

        for i = 2:m - 1
            f = f + (real(i - 1) .* tmp - 1.0e0).^2;
        end

        f = f + 1.0e0;
    elseif (strcmp(name, 'argtrig'))
        %% argtrig
        %  see j. more, 1981, testing unconstrained optimization software,;
        %  problem 26, also cuter.;
        f = 0.0e0;

        for i = 1:n
            tmp = real(n) + real(i) .* (1.0e0 - cos(y(i + 1))) - sin(y(i + 1));

            for j = 1:n
                tmp = tmp - cos(y(j + 1));
            end

            f = f + tmp .* tmp;
        end

    elseif (strcmp(name, 'arwhead'))
        %% arwhead
        f = 0.0e0;

        for i = 1:n - 1
            f = f + (y(i + 1) .* y(i + 1) + y(n + 1) .* y(n + 1)).^2 - 4.0e0 .* y(i + 1) +3.0e0;
        end
    elseif (strcmp(name, 'bdqrtic'))
        %%
        % see powell, 2003, on trust region methods for_ml unconstrained;
        % minimization without derivatives, problem 2. also see cuter.;
        % notice: the definitions in the above references are different,;
        % the difference being -4y(i)+3 v.s.(-4y(i)+3)^2. we the former;
        % one bdqrticp(see below). the original definition;
        % is in conn, gould lescrenoer, toint, 1994, performance of a;
        % multifrontal scheme for_ml partially seperable optimization, problem;
        % 61. but the paper is bitcmp available for_ml now.;
        f = 0.0e0;

        for i = 1:n - 4
            f = f + (y(i + 1).^2 + 2.0e0 .* y(i + 1 + 1).^2 + 3.0e0 .* y(i + 2 + 1).^2 + 4.0e0 .* y(i + 3 + 1).^2 + 5.0e0 .* y(n + 1).^2).^2 + (- 4.0e0 .* y(i + 1) + 3.0e0).^2;
        end
    elseif (strcmp(name, 'bdqrticp'))
        %%
        % see powell, 2003, on trust region methods for_ml unconstrained;
        % minimization without derivatives, problem 2. also see cuter.;
        % notice: the definitions in the above references are different,;
        % the difference being -4Y(i)+3 v.s. (-4Y(i)+3)^2. here is powell's;
        % version.;
        f = 0.0e0;

        for i = 1:n - 4
            f = f + (y(i + 1).^2 + 2.0e0 .* y(i + 1 + 1).^2 + 3.0e0 .* y(i + 2 + 1).^2 + 4.0e0 .* y(i + 3 + 1).^2 + 5.0e0 .* y(n + 1).^2).^2 - 4.0e0 .* y(i + 1) + 3.0e0;
        end
    elseif (strcmp(name, 'bdvalue'))
        %%
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 28, also cuter.;
        f = 0.0e0;
        h = 1.0e0 ./ real(n + 1);

        for i = 1:n
            f = f + (2.0e0 .* y(i + 1) - y(i - 1 + 1) - y(i + 1 + 1) + 0.5e0 .* h .* h .* (y(i + 1) + real(i) ./ real(n + 1) + 1.0e0).^3).^2;
        end
    elseif (strcmp(name, 'brownal'))
        %%
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 27, also cuter.;
        f = 0.0e0;
        tmp1 = 0.0e0;
        tmp2 = 1.0e0;

        for i = 1:n
            tmp1 = tmp1 + y(i + 1);
            tmp2 = tmp2 .* y(i + 1);
        end

        f = (tmp2 - 1.0e0).^2;

        for i = 1:n - 1
            f = f + (y(i + 1) + tmp1 - real(n + 1)).^2;
        end

    elseif (strcmp(name, 'broydn3d'))
        %%
        % see cuter.;
        f = 0.0e0;

        for i = 1:n
            f = f + ((3.0e0 - 2.0e0 .* y(i + 1)) .* y(i + 1) - y(i - 1 + 1) - 2 .* y(i + 1 + 1) + 1.0e0).^2;
        end

    elseif (strcmp(name, 'broydn7d'))
        %%
        % see ph. l. toint 'some numerical results using a sparse matrix;
        % updating formula in unconstrained optimization', problem 3.4, also;
        % cuter.;
        f = 0.0e0;

        for i = 1:n
            f = f + abs(y(i - 1 + 1) - y(i + 1) .* (3.0e0 - 0.5e0 .* y(i + 1)) + 2.0e0 .* y(i + 1 + 1) - 1.0e0).^(7.0e0 ./ 3.0e0);
        end

        for i = 1:fix(n ./ 2)
            f = f + (abs(y(i + 1) + y(i + fix(n ./ 2) + 1)).^(7.0e0 ./ 3.0e0));
        end

    elseif (strcmp(name, 'brybnd'))
        %%
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 31, also cuter.;
        ml = 5;
        mu = 1;
        f = 0.0e0;

        for i = 1:n
            tmp = 0.0e0;

            for j = max(1, i - ml):min(n, i + mu)

                if (j ~= i)
                    tmp = tmp + y(j + 1) .* (1.0e0 + y(j + 1));
                end

            end

            f = f + (y(i + 1) .* (2.0e0 + 5.0e00 .* y(i + 1) .* y(i + 1)) + 1.0e0 - tmp).^2;
        end

    elseif (strcmp(name, 'chainwoo'))
        %%
        % see cuter.;
        f = 1.0e0;

        for i = 1:fix(n ./ 2) - 1
            j = fix(2 .* i);
            f = f + 1.0e2 .* (y(j + 1) - y(j - 1 + 1).^2).^2 + (1.0e0 - y(j - 1 + 1)).^2 + 9.0e1 .* (y(j + 2 + 1) - y(j + 1 + 1).^2).^2 + (1.0e0 - y(j + 1 + 1)).^2 + 1.0e1 .* (y(j + 1) + y(j + 2 + 1) - 2.0e0).^2 +1.0e-1 .* (y(j + 1) - y(j + 2 + 1)).^2;
        end

    elseif (strcmp(name, 'chebquad'))
        %%
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 35, also cuter.;
        ychebq(1, :) = 1.0e0;
        ychebq(2, :) = 2.0e0 .* y([1 + 1:n + 1]) - 1.0e0;

        for i = 2:n
            ychebq(i + 1, [1:n]) = 2.0e0 .* ychebq(2, [1:n]) .* ychebq(i, [1:n]) - ychebq(i - 1, [1:n]);
        end

        f = 0.0e0;

        for i = 1:n + 1
            tmp = sum(sum(ychebq(i, :))) ./ real(n);

            if (rem(i, 2) == 1)
                tmp = tmp + 1.0e0 ./ real(i .* i - 2 .* i);
            end

            f = f + tmp .* tmp;
        end

    elseif (strcmp(name, 'chpowellb'))
        %%
        % chained version of powell badly scaled function.;
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 3.;
        f = 0.0e0;

        for i = 1:n - 1
            f = f + (1.0e4 .* y(i + 1) .* y(i + 1 + 1) - 1.0e0).^2 + (exp(-y(i + 1)) + exp(-y(i + 1 + 1)) - 1.0001e0).^2;
        end

    elseif (strcmp(name, 'chpowells'))
        %%
        % chained version of powell singular function.;
        % see j. more, 1981, testing unconstrained optimization software,;
        % problem 13, and ying-jie li, dong-hui li, truncated regularized newton;
        % method for_ml convex minimization, problem 6, and luksan, vicek,;
        % sparse and partially separable test problems for_ml unconstrained and;
        % equality constrained optimization.;
        f = 0.0e0;

        for j = 1:fix((n - 2) ./ 2)
            i = fix(2 .* j);
            f = f + (y(i - 1 + 1) + 10.0e0 .* y(i + 1)).^2 +5.0e0 .* (y(i + 1 + 1) - y(i + 2 + 1)).^2 + (y(i + 1) - 2.0e0 .* y(i + 1 + 1)).^4 +10.0e0 .* (y(i - 1 + 1) - y(i + 2 + 1)).^4;
        end

    elseif (strcmp(name, 'chnrosnb'))
        %%
        % see cuter. modified by luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.' it is another version of chained;
        % rosenbrock. compare with chrosen and rosenbrock.;
        alpha_ml = 16.0e0 .* (1.5e0 + sin(real(i))).^2;
        f = 0.0e0;

        for i = 2:n
            f = f + alpha_ml .* (y(i - 1 + 1) - y(i + 1).^2).^2 + (y(i + 1) - 1).^2;
        end

    elseif (strcmp(name, 'chrosen'))

        f = 0.0e0;

        for i = 1:n - 1
            f = f + (4.0e0) .* (y(i + 1) - y(i + 1 + 1) .* y(i + 1 + 1)) .* (y(i + 1) - y(i + 1 + 1) .* y(i + 1 + 1)) + (1.0e0 - y(i + 1 + 1)) .* (1.0e0 - y(i + 1 + 1));
        end

    elseif (strcmp(name, 'cosine'))
        %%
        % see cuter.;
        f = 0.0e0;

        for i = 1:n - 1
            f = f + cos(y(i + 1).^2 - 0.5e0 .* y(i + 1 + 1));
        end

    elseif (strcmp(name, 'cragglvy'))
        % see cuter.;
        f = 0.0e0;

        for j = 1:fix((n - 2) ./ 2)
            i = fix(2 .* j);
            f = f + (exp(y(i - 1 + 1)) - y(i + 1)).^4 ...
                + 1.0e2 .* (y(i + 1) - y(i + 1 + 1)).^6 + ...
                (tan(y(i + 1 + 1) - y(i + 2 + 1)) + y(i + 1 + 1) ...
                - y(i + 2 + 1)).^4 + y(i - 1 + 1).^8 + (y(i + 2 + 1) - 1.0e0).^2;
        end

    elseif (strcmp(name, 'cube'))
        % see cuter.;
        f = (y(1 + 1) - 1.0e0).^2;

        for i = 2:n
            f = f + 100.0e0 .* (y(i + 1) - y(i - 1 + 1).^3).^2;
        end

    elseif (strcmp(name, 'curly10') || strcmp(name, 'curly20') || strcmp(name, 'curly30'))
        %%
        % see cuter.;
        if (strcmp(name, 'curly10'))
            k = 10;
        elseif (strcmp(name, 'curly20'))
            k = 20;
        else
            k = 30;
        end

        f = 0.0e0;

        for i = 1:n
            tmp = 0.0e0;

            for j = i:min(i + k, n)
                tmp = tmp + y(j + 1);
            end

            f = f + tmp .* (tmp .* (tmp.^2 - 2.0e1) - 1.0e-1);
        end

    elseif (strcmp(name, 'dixmaane') || strcmp(name, 'dixmaanf') || strcmp(name, 'dixmaang') || strcmp(name, 'dixmaanh') || strcmp(name, 'dixmaani') || strcmp(name, 'dixmaanj') || strcmp(name, 'dixmaank') || strcmp(name, 'dixmaanl') || strcmp(name, 'dixmaanm') || strcmp(name, 'dixmaann') || strcmp(name, 'dixmaano') || strcmp(name, 'dixmaanp'))
        m = fix(fix(n ./ 3));

        if (strcmp(name, 'dixmaane') || strcmp(name, 'dixmaani') || strcmp(name, 'dixmaanm'))
            alpha_ml = 1.0e0;
            beta = 0.0e0;
            gamm = 0.125e0;
            delta = 0.125e0;
        elseif (strcmp(name, 'dixmaanf') || strcmp(name, 'dixmaanj') || strcmp(name, 'dixmaann'))
            alpha_ml = 1.0e0;
            beta = 0.0625e0;
            gamm = 0.0625e0;
            delta = 0.0625e0;
        elseif (strcmp(name, 'dixmaang') || strcmp(name, 'dixmaank') || strcmp(name, 'dixmaano'))
            alpha_ml = 1.0e0;
            beta = 0.125e0;
            gamm = 0.125e0;
            delta = 0.125e0;
        elseif (strcmp(name, 'dixmaanh') || strcmp(name, 'dixmaanl') || strcmp(name, 'dixmaanp'))
            alpha_ml = 1.0e0;
            beta = 0.26e0;
            gamm = 0.26e0;
            delta = 0.26e0;
        end

        if (strcmp(name, 'dixmaane') || strcmp(name, 'dixmaanf') || strcmp(name, 'dixmaang') || strcmp(name, 'dixmaanh'))
            p1 = 1;
            p2 = 0;
            p3 = 0;
            p4 = 1;
        elseif (strcmp(name, 'dixmaani') || strcmp(name, 'dixmaanj') || strcmp(name, 'dixmaank') || strcmp(name, 'dixmaanl'))
            p1 = 2;
            p2 = 0;
            p3 = 0;
            p4 = 2;
        elseif (strcmp(name, 'dixmaanm') || strcmp(name, 'dixmaann') || strcmp(name, 'dixmaano') || strcmp(name, 'dixmaanp'))
            p1 = 2;
            p2 = 1;
            p3 = 1;
            p4 = 2;
        end

        for i = 1:n
            f1 = f1 + (real(i) ./ real(n)).^p1 .* y(i + 1).^2;
        end

        f2 = 0.0e0;

        for i = 1:n - 1
            f2 = f2 + (real(i) ./ real(n)).^p2 .* y(i + 1).^2 .* (y(i + 1 + 1) + y(i + 1 + 1).^2).^2;
        end

        f3 = 0.0e0;

        for i = 1:2 * m
            f3 = f3 + (real(i) ./ real(n)).^p3 .* y(i + 1).^2 .* y(i + m + 1).^4;
        end

        f4 = 0.0e0;

        for i = 1:m
            f4 = f4 + (real(i) ./ real(n)).^p4 .* y(i + 1) .* y(i + 2 .* m + 1);
        end

        f = 1.0e0;
        f = f + alpha_ml .* f1 + beta .* f3 + gamm .* f3 + delta .* f4;
    elseif (strcmp(name, 'dqrtic'))
        %%
        % c see cuter.;
        f = 0.0e0;

        for i = 1:n
            f = f + (y(i + 1) - real(i)).^4;
        end

    elseif (strcmp(name, 'edensch'))
        %%
        % see cuter.;
        f = 16.0e0;

        for i = 1:n - 1
            f = f + (y(i + 1) - 2.0e0).^4 + (y(i + 1) .* y(i + 1 + 1) - 2.0e0 .* y(i + 1 + 1)).^2 + (y(i + 1 + 1) + 1.0e0).^2;
        end

    elseif (strcmp(name, 'eg2'))
        %%
        % see cute;
        f = 0.0e0;

        for i = 1:n - 1
            f = f + sin(y(1 + 1) + y(i + 1).^2 - 1.0e0);
        end

        f = f + 0.5e0 .* sin(y(n + 1).^2);
    elseif (strcmp(name, 'engval1'))
        %%
        % see cuter. compare with arwhead.;
        f = 0.0e0;

        for i = 1:n - 1
            f = f + (y(i + 1).^2 + y(i + 1 + 1).^2).^2 - 4.0e0 .* y(i + 1) + 3.0e0;
        end

    elseif (strcmp(name, 'errinros'))
        %%
        % see cuter. modified by luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.' it is another version of chained;
        % rosenbrock. compare with chnrosnb, chrosen, and rosenbrock.;
        f = 0.0e0;

        for i = 2:n
            alpha = 16 * (1.5 + sin(i))^2;
            f = f + (y(i) - alpha .* y(i + 1)^2)^2 + (y(i + 1) - 1)^2;
        end

    elseif (strcmp(name, 'expsum'))
        %%
        % see http:[,www.opt.uni]-duesseldorf.de./~jarre./dot./f_cx.m;
        alpha_ml = 4.0e0;
        beta = 2.0e0 .* sqrt(alpha_ml);
        f = 0.0e0;

        for i = 1:n
            f = f + real(i.^2) .* (exp(real(i) .* sum(sum(y([1 + 1:i + 1])))) + alpha_ml .* exp(-real(i) .* sum(sum(y([1 + 1:i + 1])))) - beta);
        end

    elseif (strcmp(name, 'extrosnb'))
        % see cuter and luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.' it is another version of chained;
        % rosenbrock. compare with chnrosnb, chrosen, errinros, and rosenbrock.;
        f = (y(1 + 1) - 1.0e0).^2;

        for i = 2:n
            f = f + 1.0e2 .* (y(i + 1) - y(i - 1 + 1).^2).^2;
        end

    elseif (strcmp(name, 'exttet'))
        % see n. andrei, 2008, an unconstrained optimization test functions;
        % collection, extended three term exponentials.;
        f = 0.0e0;

        for j = 1:fix(n ./ 2)
            i = fix(2 .* j);
            f = f + exp(y(i - 1 + 1) + 3.0e0 .* y(i + 1) - 0.1e0) + exp(y(i - 1 + 1) - 3.0e0 .* y(i + 1) - 0.1e0) + exp(-y(i - 1 + 1) - 0.1e0);
        end

    elseif (strcmp(name, 'firose'))
        % five-diagonal rosenbrock. the jacobian of the corresponding nonlinear;
        % equations is five-diagonal. see example 5.3 of 'the secant/finite;
        % difference algorithms for_ml solving sparse nonlinear systems of equations';
        % by guangye li(siam journal on numerical analysis, 1988).;
        f = 0.0e0;

        if (n >= 4)
            f = (4.0e0 .* (y(1 + 1) - y(2 + 1).^2) + y(2 + 1) - y(3 + 1).^2).^2 + (8.0e0 .* y(2 + 1) .* (y(2 + 1).^2 - y(1 + 1)) - 2.0e0 .* (1.0e0 - y(2 + 1)) + 4.0e0 .* (y(2 + 1) - y(3 + 1).^2) + y(3 + 1) - y(4 + 1).^2).^2;

            for i = 3:n - 2
                f = f + (8.0e0 .* y(i + 1) .* (y(i + 1).^2 - y(i - 1 + 1)) - 2.0e0 .* (1.0e0 - y(i + 1)) + 4.0e0 .* (y(i + 1) - y(i + 1 + 1).^2) + y(i - 1 + 1).^2 - y(i - 2 + 1) + y(i + 1 + 1) - y(i + 2 + 1).^2).^2;
            end

            f = f + (8.0e0 .* y(n - 1 + 1) .* (y(n - 1 + 1).^2 - y(n - 2 + 1)) - 2.0e0 .* (1.0e0 - y(n - 1 + 1)) + 4.0e0 .* (y(n - 1 + 1) - y(n + 1).^2) + y(n - 2 + 1).^2 - y(n - 3 + 1)).^2;
            f = f + (8.0e0 .* y(n + 1) .* (y(n + 1).^2 - y(n - 1 + 1)) - 2.0e0 .* (1.0e0 - y(n + 1)) + y(n - 1 + 1).^2 - y(n - 2 + 1)).^2;
        end

    elseif (strcmp(name, 'fletcbv2'))
        %%
        % see cuter.
        h = 1.0e0 ./ real(n + 1);
        tmp1 = y(1 + 1).^2;

        for i = 1:n - 1
            tmp1 = tmp1 + (y(i + 1) - y(i + 1 + 1)).^2;
        end

        tmp1 = tmp1 + y(n + 1).^2;
        tmp2 = 0.0e0;

        for i = 1:n
            tmp2 = tmp2 + 2.0e0 .* y(i + 1) + cos(y(i + 1));
        end

        f = 0.5e0 .* tmp1 - h.^2 .* tmp2 -y(n + 1);
    elseif (strcmp(name, 'fletcbv3'))
        %%
        % see cuter. modified by luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.';
        p = 1.0e-8;
        h = 1.0e0 ./ real(n + 1);
        kappa = 1.0e0;
        tmp1 = y(1 + 1).^2;

        for i = 1:n - 1
            tmp1 = tmp1 + (y(i + 1) - y(i + 1 + 1)).^2;
        end

        tmp1 = tmp1 + y(n + 1).^2;
        tmp2 = 0.0e0;

        for i = 1:n
            tmp2 = tmp2 + 1.0e2 .* sin(y(i + 1) .* 1.0e-2) .* (1.0e0 + 2.0e0 ./ h.^2) + 1.0e0 ./ h.^2 .* kappa .* cos(y(i + 1));
        end

        f = 0.5e0 .* p .* tmp1 - p .* tmp2;
    elseif (strcmp(name, 'fletchcr'))
        %%
        % see cuter.;
        f = 0.0e0;

        for i = 1:n - 1
            f = f + (y(i + 1 + 1) - y(i + 1) + 1.0e0 - y(i + 1).^2).^2;
        end

        f = 1.0e2 .* f;
    elseif (strcmp(name, 'fminsrf2'))
        %%
        % see cuter;
        f = 0.0e0;
        k = floor(sqrt(real(n)));
        hk = floor(real(k) ./ 2.0e0);

        if (k <= 1)
            info = 2;
            f = 0;
        else

            for i = 1:k

                for j = 1:k
                    ysrf(i, j) = x((i - 1) .* k + j);
                end

            end

            for i = 1:k - 1

                for j = 1:k - 1
                    f = f + sqrt(1.0e0 + 0.5e0 .* real(k - 1).^2 .* (((ysrf(i, j) - ysrf(i + 1, j + 1)).^2 + (ysrf(i + 1, j) - ysrf(i, j + 1)).^2)));
                end

            end

            f = f ./ real(k - 1).^2 + ysrf(hk, hk).^2 ./ real(n);
        end

    elseif (strcmp(name, 'freuroth'))
        %%
        tmp1 = 0.0e0;
        tmp2 = 0.0e0;

        for i = 1:n - 1
            tmp1 = tmp1 + ((5.0e0 - y(i + 1 + 1)) .* y(i + 1 + 1).^2 + y(i + 1) - 2.0e0 .* y(i + 1 + 1) - 13.0e0).^2;
            tmp2 = tmp2 + ((1.0e0 + y(i + 1 + 1)) .* y(i + 1 + 1).^2 + y(i + 1) - 14.0e0 .* y(i + 1 + 1) - 29.0e0).^2;
        end

        f = tmp1 + tmp2;
    elseif (strcmp(name, 'genbrown'))
        %%
        % generalized brown function 1. see ying-jie li, dong-hui li, truncated;
        % regularized newton method for_ml convex minimization, problem 6,;
        % and luksan, vicek, sparse and partially separable test problems for_ml;
        % unconstrained and equality constrained optimization.;
        f = 0.0e0;

        for i = 2:n
            f = f + (y(i - 1 + 1) - 3.0e0).^2 + (y(i - 1 + 1) - y(i + 1)).^2 ...
                + exp(20.0e0 .* (y(i - 1 + 1) - y(i + 1)));
        end

    elseif (strcmp(name, 'genhumps'))
        %%
        % see cuter.;
        zeta = 2.0e0;
        f = 0.0e0;

        for i = 1:n - 1
            f = f + sin(zeta .* y(i + 1)).^2 .* sin(zeta .* y(i + 1 + 1)).^2 + 0.05e0 .* (y(i + 1).^2 + y(i + 1 + 1).^2);
        end

    elseif (strcmp(name, 'genrose'))
        %%
        % see cuter.;
        % compare with chrosen and rosenbrock.;
        % the starting point is(fix(i./(n+1))).;
        f = 1.0e0;

        for i = 2:n
            f = f + 1.0e2 .* (y(i + 1) - y(i - 1 + 1).^2).^2 + (y(i + 1) - 1.0e0).^2;
        end

    elseif (strcmp(name, 'indef'))
        %%
        % see cuter. modified by luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.';
        tmp1 = 0.0e0;

        for i = 1:n
            tmp1 = tmp1 + 1.0e2 .* sin(1.0e-2 .* y(i + 1));
        end

        tmp2 = 0.0e0;

        for i = 2:n - 1
            tmp2 = tmp2 + cos(2.0e0 .* y(i + 1) - y(n + 1) - y(1 + 1));
        end

        f = tmp1 + 0.5e0 .* tmp2;
    elseif (strcmp(name, 'integreq'))
        %%
        % see j. more, 1981, testing unconstrained optimization software,;
        %  problem 29 and also cuter.;
        h = 1.0e0 ./ real(n + 1);

        for i = 1:n
            t(i) = real(i) ./ real(n + 1);
        end

        f = 0.0e0;

        for i = 1:n
            tmp1 = 0.0e0;

            for j = 1:i
                tmp1 = tmp1 + t(j) .* (y(j + 1) + t(j) + 1.0e0).^3;
            end

            tmp2 = 0.0e0;

            for j = i + 1:n
                tmp2 = tmp2 + (1.0e0 - t(j)) .* (y(j + 1) + t(j) + 1.0e0).^3;
            end

            f = f + (y(i + 1) + 0.5e0 .* h .* ((1 - t(i)) .* tmp1 + t(i) .* tmp2)).^2;
        end

    elseif (strcmp(name, 'liarwhd'))
        %%
        f = 0.0e0;

        for i = 1:n
            f = f + 4.0e0 .* (y(i + 1).^2 - y(1 + 1)).^2 + (y(i + 1) - 1.0e0).^2;
        end

    elseif (strcmp(name, 'lilifun3'))
        %%
        % see ying-jie li, dong-hui li, truncated;
        % regularized newton method for_ml convex minimization, problem 3.;
        f = 0.0e0;

        for i = 2:n
            f = f + exp(y(i + 1) - y(i - 1 + 1)).^2 + (y(i + 1) - y(i - 1 + 1)).^2 ...
                + 2.0e0 .* y(i + 1).^4 + 4.0e0 .* y(i - 1 + 1).^4;
        end

    elseif (strcmp(name, 'lilifun4'))
        % see ying-jie li, dong-hui li, truncated;
        % regularized newton method for_ml convex minimization, problem 4.;
        f = 0.0e0;

        for i = 2:n
            f = f + 0.5e0 .* (y(i + 1) - y(i - 1 + 1)).^2 + sin(y(i + 1) - y(i - 1 + 1)) ...
                + 2.0e0 * (2.0e0 .* y(i + 1) + 3.0e0 .* y(i - 1 + 1) - 15.0e0).^4;
        end

    elseif (strcmp(name, 'morebv') || strcmp(name, 'morebvl'))
        %%
        % see cuter. morebvl uses the start point suggested by luksan, matonoha,;
        % vlcek, 'modified cute probelems for;
        % unconstrained optimzation.';
        h = 1.0e0 ./ real(n + 1);
        f = 0.0e0;

        for i = 1:n
            f = f + (2.0e0 .* y(i + 1) - y(i - 1 + 1) - y(i + 1 + 1) + 0.5e0 .* h .* h .* (y(i + 1) + real(i) .* h + 1.0e0).^3).^2;
        end

    elseif (strcmp(name, 'ncb20'))
        %%
        % see cuter. corrected by luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.';
        f = 0.0e0;

        if (n >= 20)
            tmp1 = 0.0e0;

            for i = 1:n - 30
                tmp = 0.0e0;

                for j = 1:20
                    tmp = tmp + y(i + j - 1 + 1) ./ (1 + y(i + j - 1 + 1).^2);
                end

                tmp1 = tmp1 + 10.0e0 ./ real(i) .* tmp.^2;
                tmp = 0.0e0;

                for j = 1:20
                    tmp = tmp + y(i + j - 1 + 1);
                end

                tmp1 = tmp1 - 0.2e0 .* tmp;
            end

            tmp2 = 0.0e0;

            for i = 1:n - 10
                tmp2 = tmp2 + y(i + 1).^4 + 2.0e0;
            end

            tmp3 = 0.0e0;

            for i = 1:10
                tmp3 = tmp3 + y(i + 1) .* y(i + 10 + 1) .* y(i + n - 10 + 1) +2.0e0 .* y(i + n - 10 + 1).^2;
            end

            f = 2.0e0 + tmp1 + tmp2 + 1.0e-4 .* tmp3;
        end

    elseif (strcmp(name, 'ncb20b'))
        %%
        % see cuter. corrected by luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.';
        f = 0.0e0;

        if (n >= 20)
            tmp1 = 0.0e0;

            for i = 1:n - 19
                tmp = 0.0e0;

                for j = 1:20
                    tmp = tmp + y(i + j - 1 + 1) ./ (1 + y(i + j - 1 + 1).^2);
                end

                tmp1 = tmp1 + 10.0e0 ./ real(i) .* tmp.^2;
                tmp = 0.0e0;

                for j = 1:20
                    tmp = tmp + y(i + j - 1 + 1);
                end

                tmp1 = tmp1 - 0.2e0 .* tmp;
            end

            tmp2 = 0.0e0;

            for i = 1:n
                tmp2 = tmp2 + (1.0e2 .* y(i + 1).^4 + 2.0e0);
            end

            f = tmp1 + tmp2;
        end

    elseif (strcmp(name, 'noncvxu2'))
        %%
        % see cuter.;
        f = 0.0e0;

        for i = 1:n
            j = fix(rem(3 .* i - 2, n));
            k = fix(rem(7 .* i - 3, n));
            tmp = y(i + 1) + y(j + 1 + 1) + y(k + 1 + 1);
            f = f + tmp.^2 + 4.0e0 .* cos(tmp);
        end

    elseif (strcmp(name, 'noncvxun'))
        % see cuter.;
        f = 0.0e0;

        for i = 1:n
            j = fix(rem(2 .* i - 1, n));
            k = fix(rem(3 .* i - 1, n));
            tmp = y(i + 1) + y(j + 1 + 1) + y(k + 1 + 1);
            f = f + tmp.^2 + 4.0e0 .* cos(tmp);
        end

    elseif (strcmp(name, 'nondia'))
        % see cuter.;
        f = (y(1 + 1) - 1.0e0).^2;

        for i = 2:n
            f = f + (1.0e2 .* y(1 + 1) - y(i - 1 + 1).^2).^2;
        end

    elseif (strcmp(name, 'nondquar'))
        f = (y(1 + 1) - y(2 + 1)).^2 + (y(n - 1 + 1) - y(n + 1)).^2;

        for i = 1:n - 2
            f = f + (y(i + 1) + y(i + 1 + 1) + y(n + 1)).^4;
        end

    elseif (strcmp(name, 'penalty1'))
        dp1 = 0.0e0;
        dp2 = 0.0e0;

        for i = 1:n
            dp1 = dp1 + (y(i + 1) - 1.0e0).^2;
            dp2 = dp2 + y(i + 1).^2;
        end

        f = 1.0e-5 .* dp1 + (0.25e0 - dp2).^2;
    elseif (strcmp(name, 'penalty2'))
        %%
        f = 0.0e0;
        tmp = 0.0e0;

        for i = 2:n
            f = f + (exp(y(i - 1 + 1) .* 0.1e0) + exp(y(i + 1) .* 0.1e0) - exp(real(i - 1) .* 0.1e0) - exp(real(i) .* 0.1e0)).^2 + (exp(y(i + 1) .* 0.1e0) -exp(-0.1e0)).^2;
        end

        for i = 1:n
            tmp = tmp + real(n - i + 1) .* y(i + 1) .* y(i + 1);
        end

        f = f + (1.0e0 - tmp) .* (1.0e0 - tmp) + (y(1 + 1) - 0.2e0) .* (y(1 + 1) - 0.2e0);
    elseif (strcmp(name, 'penalty3') || strcmp(name, 'penalty3p'))
        %%
        f = 0.0e0;
        r = 0.0e0;
        s = 0.0e0;
        pd = 0.0e0;

        for i = 1:n
            pd = pd + (y(i + 1).^2 - real(n));
        end

        if (strcmp(name, 'penalty3'))
            alpha_ml = 1.0e0;
        else
            % Suggested by Powell, 'The NEWUOA Software
            alpha_ml = 1.0e-3;
            %           for_ml unconstrained optimization without derivatives', 2004,;
            %            section 8.;
        end

        for i = 1:n - 2
            r = r + (y(i + 1) + 2.0e0 .* y(i + 1 + 1) + 1.0e1 .* y(i + 2 + 1) - 1.0e0).^2;
            s = s + (2.0e0 .* y(i + 1) + y(i + 1 + 1) - 3.0e0).^2;
        end

        for i = 1:fix(n ./ 2)
            f = f + (y(i + 1) - 1.0e0).^2;
        end

        f = f + pd.^2 + alpha_ml .* (1.0e0 + r .* exp(y(n + 1)) + s .* exp(y(n - 1 + 1)) +r .* s);
    elseif (strcmp(name, 'powellsg'))
        %%
        f = 0.0e0;

        for i = 1:fix(n ./ 4)
            j = fix(4 .* (i - 1) + 1);
            f = f + (y(j + 1) + 1.0e1 .* y(j + 1 + 1)).^2 +5.0e0 .* (y(j + 2 + 1) - y(j + 3 + 1)).^2 + (y(j + 1 + 1) - 2.0e0 .* y(j + 2 + 1)).^4 + 1.0e1 .* (y(j + 1) - y(j + 3 + 1)).^4;
        end

    elseif (strcmp(name, 'power'))
        %%
        % see cuter.;
        f = 0.0e0;

        for i = 1:n
            f = f + (real(i) .* y(i + 1)).^2;
        end

    elseif (strcmp(name, 'rosenbrock'))
        %%
        % in claasical rosenbrock function, alpha_ml = 100.0e0;
        % when alpha_ml = 4.0e0, this function is essentially the same as chrosen,;
        % except the order of the variables.;
        alpha_ml = 100.0e0;
        f = 0.0e0;

        for i = 1:n - 1
            f = f + (1.0e0 - y(i + 1)) .* (1.0e0 - y(i + 1)) + alpha_ml .* (y(i + 1 + 1) - y(i + 1) .* y(i + 1)) .* (y(i + 1 + 1) - y(i + 1) .* y(i + 1));
        end

    elseif (strcmp(name, 'sbrybnd') || strcmp(name, 'sbrybndl'))
        %%
        % see cuter. and see luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.' the scaling parameter of luksan is 6;
        % instead of 12. there is a typo in luksan: the square is missed, which;
        % makes the function unbounded from below.;
        ml = 5;
        mu = 1;

        if (strcmp(name, 'sbrybnd'))
            scaling = 12.0e0;
        else
            scaling = 6.0e0;
        end

        f = 0.0e0;

        if (n >= 2)

            for i = 1:n
                ytmp(i) = exp(scaling .* (real(i - 1) ./ real(n - 1))) .* y(i + 1);
            end

            for i = 1:n
                tmp = 0.0e0;

                for j = max(1, i - ml):min(n, i + mu)

                    if (j ~= i)
                        tmp = tmp + ytmp(j) .* (1.0e0 + ytmp(j));
                    end

                end

                f = f + (ytmp(i) .* (2.0e0 + 5.0e00 .* ytmp(i) .* ytmp(i)) + 1.0e0 - tmp).^2;
            end

        end

    elseif (strcmp(name, 'schmvett'))
        % see cuter.;
        f = 0.0e0;

        for i = 1:n -2
            f = f - 1.0e0 ./ (1.0e0 + (y(i + 1) - y(i + 1 + 1)).^2) - sin(0.5e0 .* (pi .* y(i + 1 + 1) + y(i + 2 + 1))) - exp(- ((y(i + 1) + y(i + 2 + 1)) ./ y(i + 1 + 1) - 2.0e0).^2);
        end

    elseif (strcmp(name, 'scosine') || strcmp(name, 'scosinel'))
        %%
        % see cuter. and see luksan, matonoha, vlcek, 'modified cute probelems for;
        % unconstrained optimzation.' the scaling parameter of luksan is 6;
        % instead of 12.;
        if (strcmp(name, 'scosine'))
            scaling = 12.0e0;
        else
            scaling = 6.0e0;
        end

        f = 0.0e0;

        if (n >= 2)

            for i = 1:n
                ytmp(i) = exp(scaling .* (real(i - 1) ./ real(n - 1))) .* y(i + 1);
            end

            for i = 1:n - 1
                f = f + cos(ytmp(i).^2 - 0.5e0 .* ytmp(i + 1));
            end

        end

    elseif (strcmp(name, 'serose'))
        %%
        % seven-diagonal rosenbrock. the jacobian of the corresponding nonlinear;
        % equations is 7-diagonal. see example 5.4 of 'the secant/finite;
        % difference algorithms for_ml solving sparse nonlinear systems of equations';
        % by guangye li(siam journal on numerical analysis, 1988).;
        % the formulation in the paper seems bitcmp correct. the cases with j=3 and;
        % j=n-2 should be defined seperately as well.;
        f = 0.0e0;

        if (n >= 6)
            f = (4.0e0 .* (y(1 + 1) - y(2 + 1).^2) + y(2 + 1) - y(3 + 1).^2 + y(3 + 1) - y(4 + 1).^2).^2;
            f = f + (8.0e0 .* y(2 + 1) .* (y(2 + 1).^2 - y(1 + 1)) - 2.0e0 .* (1.0e0 - y(2 + 1)) + 4.0e0 .* (y(2 + 1) - y(3 + 1).^2) + y(3 + 1) - y(4 + 1).^2 + y(4 + 1) - y(5 + 1).^2).^2;
            f = f + (8.0e0 .* y(3 + 1) .* (y(3 + 1).^2 - y(2 + 1)) - 2.0e0 .* (1.0e0 - y(3 + 1)) + 4.0e0 .* (y(3 + 1) - y(4 + 1).^2) + y(2 + 1).^2 - y(1 + 1) + y(4 + 1) - y(5 + 1).^2 + y(5 + 1) - y(6 + 1).^2).^2;

            for i = 4:n - 3
                f = f + (8.0e0 .* y(i + 1) .* (y(i + 1).^2 - y(i - 1 + 1)) - 2.0e0 .* (1.0e0 - y(i + 1)) + 4.0e0 .* (y(i + 1) - y(i + 1 + 1).^2) + y(i - 1 + 1).^2 - y(i - 2 + 1) + y(i + 1 + 1) - y(i + 2 + 1).^2 + y(i - 2 + 1).^2 - y(i - 3 + 1) + y(i + 2 + 1) - y(i + 3 + 1).^2).^2;
            end

            f = f + (8.0e0 .* y(n - 2 + 1) .* (y(n - 2 + 1).^2 - y(n - 3 + 1)) - 2.0e0 .* (1.0e0 - y(n - 2 + 1)) + 4.0e0 .* (y(n - 2 + 1) - y(n - 1 + 1).^2) + y(n - 3 + 1).^2 - y(n - 4 + 1) + y(n - 1 + 1) - y(n + 1).^2 + y(n - 4 + 1).^2 - y(n - 5 + 1)).^2;
            f = f + (8.0e0 .* y(n - 1 + 1) .* (y(n - 1 + 1) - y(n - 2 + 1)) - 2.0e0 .* (1.0e0 - y(n - 1 + 1)) + 4.0e0 .* (y(n - 1 + 1) - y(n + 1).^2) + y(n - 2 + 1).^2 - y(n - 3 + 1) + y(n - 3 + 1).^2 - y(n - 4 + 1)).^2;
            f = f + (8.0e0 .* y(n + 1) .* (y(n + 1).^2 - y(n - 1 + 1)) - 2.0e0 .* (1.0e0 - y(n + 1)) + y(n - 1 + 1).^2 - y(n - 2 + 1) + y(n - 2 + 1).^2 - y(n - 3 + 1)).^2;
        end

    elseif (strcmp(name, 'sinquad'))
        % see cuter.;
        f = (y(1 + 1) - 1.0e0).^4;

        for i = 2:n - 1
            f = f + (sin(y(i + 1) - y(n + 1)) - y(1 + 1).^2 + y(i + 1).^2).^2;
        end

        f = f + (y(n + 1).^2 - y(1 + 1).^2).^2;
    elseif (strcmp(name, 'sphrpts'))
        f = 0.0e0;

        if (rem(n, 2) ~= 0)
            info = 2;
        else

            for i = 2:fix(n ./ 2)

                for j = 1:i - 1
                    f = f + 1.0e0 ./ ((cos(y(2 .* i - 1 + 1)) .* cos(y(2 .* i + 1)) - cos(y(2 .* j - 1 + 1)) .* cos(y(2 .* j + 1))).^2 + (sin(y(2 .* i - 1 + 1)) .* cos(y(2 .* i + 1)) - sin(y(2 .* j - 1 + 1)) .* cos(y(2 .* j + 1))).^2 + (sin(y(2 .* i + 1)) - sin(y(2 .* j + 1))).^2);
                end

            end

        end

    elseif (strcmp(name, 'sparsine'))
        f = 0.0e0;

        for i = 1:n
            j1 = fix(rem(2 .* i - 1, n) + 1);
            j2 = fix(rem(3 .* i - 1, n) + 1);
            j3 = fix(rem(5 .* i - 1, n) + 1);
            j4 = fix(rem(7 .* i - 1, n) + 1);
            j5 = fix(rem(11 .* i - 1, n) + 1);
            f = f + real(i) .* (sin(y(i + 1)) + sin(y(j1 + 1)) + sin(y(j2 + 1)) + sin(y(j3 + 1)) + sin(y(j4 + 1)) + sin(y(j5 + 1))).^2;
        end

        f = 0.5e0 .* f;
    elseif (strcmp(name, 'sparsqur'))
        %%
        f = 0.0e0;

        for i = 1:n
            j1 = fix(rem(2 .* i - 1, n) + 1);
            j2 = fix(rem(3 .* i - 1, n) + 1);
            j3 = fix(rem(5 .* i - 1, n) + 1);
            j4 = fix(rem(7 .* i - 1, n) + 1);
            j5 = fix(rem(11 .* i - 1, n) + 1);
            f = f + real(i) .* (y(i + 1).^2 + y(j1 + 1).^2 + y(j2 + 1).^2 + y(j3 + 1).^2 + y(j4 + 1).^2 + y(j5 + 1).^2).^2;
        end

        f = 0.125e0 .* f;
    elseif (strcmp(name, 'spmsrtls'))
        %%
        for i = 1:n
            par(i) = sin(real(i).^2);
        end

        m = fix(fix((n + 2) ./ 3));
        f = 0.0e0;

        for i = 1:m
            f1 = 0.0e0;
            f2 = 0.0e0;
            f3 = 0.0e0;
            f4 = 0.0e0;
            f5 = 0.0e0;
            j = fix(3 .* (i - 1) + 1);

            if (i > 2)
                f1 = (y(j - 4 + 1) .* y(j - 1 + 1) - par(j - 4) .* par(j - 1)).^2;
            end

            f3 = (y(j + 1).^2 - par(j).^2).^2;

            if (i > 1)
                f2 = (y(j - 3 + 1) .* y(j - 1 + 1) + y(j - 1 + 1) .* y(j + 1) - par(j - 3) .* par(j - 1) - par(j - 1) .* par(j)).^2;
                f3 = f3 + (y(j - 2 + 1) .* y(j - 1 + 1) - par(j - 2) .* par(j - 1)).^2;
            end

            if (i < m)
                f3 = f3 + (y(j + 2 + 1) .* y(j + 1 + 1) - par(j + 2) .* par(j + 1)).^2;
                f4 = (y(j + 3 + 1) .* y(j + 1 + 1) + y(j + 1 + 1) .* y(j + 1) - par(j + 3) .* par(j + 1) - par(j + 1) .* par(j)).^2;
            end

            if (i < m - 1)
                f5 = (y(j + 4 + 1) .* y(j + 1 + 1) - par(j + 4) .* par(j + 1)).^2;
            end

            f = f + f1 + f2 + f3 + f4 + f5;
        end

    elseif (strcmp(name, 'srosenbr'))
        %%
        f = 0.0e0;

        for i = 1:fix(n ./ 2)
            f = f + 1.0e2 .* (y(2 .* i + 1) - y(2 .* i - 1 + 1).^2).^2 ...
                + (y(2 .* i - 1 + 1) - 1.0e0).^2;
        end

    elseif (strcmp(name, 'stmod'))
        %%
        % modified stybilsky-tang-function.;
        % see http:[,www.opt.uni]-duesseldorf.de./~jarre./dot./f_stmod.m;
        f = 0.0e0;

        for i = 1:n
            MST((1:i), i) = real(i);
            MST(i, i) = MST(i, i) + 1.0e0;
        end

        y([1 + 1:n + 1]) = (MST * y([1 + 1:n + 1]));

        for i = 1:n
            f = f + y(i + 1).^4 - 1.6e1 .* y(i + 1).^2 + 3.5e1 .* y(i + 1);
        end

        f = f + real(n) .* 1.7119854821721323e2;
    elseif (strcmp(name, 'tointgss'))
        %%
        f = 0.0e0;

        for i = 1:n - 2
            f = f + (1.0e1 ./ real(n + 2) + y(i + 2 + 1).^2) .* (2.0e0 - exp(- (y(i + 1) - y(i + 1 + 1)).^2 ./ (1.0e-1 + y(i + 2 + 1).^2)));
        end

    elseif (strcmp(name, 'tointtrig'))
        %%
        % see ph. l. toint, 'some numerical results using a sparse matrix;
        % updating formula in unconstrained optimization', example 3.6.;
        f = 0.0e0;

        for i = 2:n
            betai = 1.0e0 + 1.0e-1 .* real(i);

            for j = 1:i - 1
                alpha_ml = real(5 .* (1 + rem(i, 5) + rem(j, 5)));
                betaj = 1.0e0 + 1.0e-1 .* real(j);
                c = real(i + j) .* 1.0e-1;
                f = f + alpha_ml .* sin(betai .* y(i + 1) + betaj .* y(j + 1) + c);
            end

        end

    elseif (strcmp(name, 'tquartic'))
        %%
        f = (y(1 + 1) - 1.0e0).^2;

        for i = 1:n - 1
            f = f + (y(1 + 1).^2 - y(i + 1 + 1).^2).^2;
        end

    elseif (strcmp(name, 'trigsabs'))
        %%
        % see powell, 'The NEWUOA Software for unconstrained optimization without derivatives', 2004, section 8.;
        % the random numbers are replaced by 'irregular' numbers generated with trigonometrics.;
        for i = 1:2 * n

            for j = 1:n
                MaS(i, j) = 1.0e2 .* sin(real(i .* j .* n .* 163).^(2.5e-1));
                MaC(i, j) = 1.0e2 .* cos(real(i + j + n + 42) ./ 4.0e0);
            end

        end

        for j = 1:n
            sinx(j) = sin(y(j + 1));
            cosx(j) = cos(y(j + 1));
            so(j) = sin(pi .* cos(real(n + j + 163) ./ 3.0e0));
            co(j) = cos(pi .* cos(real(n + j + 163) ./ 3.0e0));
        end

        d = MaS * (sinx - so) + MaC * (cosx - co);
        f = sum(sum(abs(d)));
    elseif (strcmp(name, 'trigssqs'))
        %%
        % see powell, 'The NEWUOA Software for unconstrained optimization without derivatives', 2004, section 8.;
        % the random numbers are replaced by 'irregular' numbers generated with trigonometrics.;
        for i = 1:2 * n

            for j = 1:n
                MaS(i, j) = 1.0e2 .* sin(real(i .* j .* n .* 163).^(2.5e-1));
                MaC(i, j) = 1.0e2 .* cos(real(i + j + n + 42) ./ 4.0e0);
            end

        end

        for j = 1:n
            theta = (1.0e0 + sin(real(n .* j .* 42).^(1.0e0 ./ 3.0e0))) ./ 2.0e0;
            theta = 1.0e1.^(-theta);
            sinx(j) = sin(theta .* y(j + 1));
            cosx(j) = cos(theta .* y(j + 1));
            so(j) = sin(pi .* cos(real(n + j + 163) ./ 3.0e0));
            co(j) = cos(pi .* cos(real(n + j + 163) ./ 3.0e0));
        end

        d = MaS * (sinx - so) + MaC * (cosx - co);
        f = dot(d, d);
    elseif (strcmp(name, 'trirose1'))
        %%
        % tri-diagonal rosenbrock. the jacobian of the corresponding nonlinear;
        % equations is 3-diagonal. see example 5.1 of 'the secant/finite;
        % difference algorithms for_ml solving sparse nonlinear systems of equations';
        % by guangye li(siam journal on numerical analysis, 1988).;
        f = 0.0e0;

        if (n >= 2)
            f = 64.0e0 .* (y(1 + 1) - y(2 + 1).^2).^2;

            for i = 2:n - 1
                f = f + (16.0e0 .* y(i + 1) .* (y(i + 1).^2 - y(i - 1 + 1)) - 2.0e0 .* (1.0e0 - y(i + 1)) + 8.0e0 .* (y(i + 1) - y(i + 1 + 1).^2)).^2;
            end

            f = f + (16.0e0 .* y(n + 1) .* (y(n + 1).^2 - y(n - 1 + 1)) - 2.0e0 .* (1.0e0 - y(n + 1))).^2;
        end

    elseif (strcmp(name, 'trirose2'))
        %%
        % tri-diagonal rosenbrock. the jacobian of the corresponding nonlinear;
        % equations is 3-diagonal. see example 5.2 of 'the secant/finite;
        % difference algorithms for_ml solving sparse nonlinear systems of equations';
        % by guangye li(siam journal on numerical analysis, 1988).;
        f = 0.0e0;

        if (n >= 2)
            f = 16.0e0 .* (y(1 + 1) - y(2 + 1).^2).^2;

            for i = 2:n - 1
                f = f + (8.0e0 .* y(i + 1) .* (y(i + 1).^2 - y(i - 1 + 1)) - 2.0e0 .* (1.0e0 - y(i + 1)) + 4.0e0 .* (y(i + 1) - y(i + 1 + 1).^2)).^2;
            end

            f = f + (8.0e0 .* y(n + 1) .* (y(n + 1).^2 - y(n - 1 + 1)) - 2.0e0 .* (1.0e0 - y(n + 1))).^2;
        end

    elseif (strcmp(name, 'vardim'))
        f = 0.0e0;

        for i = 1:n
            f = f + real(i) .* (y(i + 1) - 1.0e0);
        end

        f = f .* f;
        f = f + f .* f;

        for i = 1:n
            f = f + (y(i + 1) - 1.0e0) .* (y(i + 1) - 1.0e0);
        end

    elseif (strcmp(name, 'woods'))
        f = 0.0e0;

        for i = 1:fix(n ./ 4)
            j = fix(4 .* i);
            f = f + 1.0e2 .* (y(j - 2 + 1) - y(j - 3 + 1).^2).^2 + (1.0e0 - y(j - 3 + 1)).^2 + 9.0e1 .* (y(j + 1) - y(j - 1 + 1).^2).^2 + (1.0e0 - y(j - 1 + 1)).^2 +1.0e1 .* (y(j - 2 + 1) + y(j + 1) - 2.0e0).^2 + 1.0e-1 .* (y(j - 2 + 1) - y(j + 1)).^2;
        end

    else
        % info = 1 means unknown function name
        info = 1;
    end

    fun = f;

end
