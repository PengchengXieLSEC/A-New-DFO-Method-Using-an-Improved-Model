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
function [b, c, itfun] = trust_sub_rfzero(x, itbnd, eigval, alpha, delta, tol)

% start the iteration counter
    itfun = 0;

    % find the first three points, a, b, and c and their values
    if (x ~= 0.0)
        dx = abs(x) / 2;
    else
        dx = 0.5;
    end

    a = x;
    c = a;
    fa = trust_sub_secular_eqn(a, eigval, alpha, delta);
    itfun = itfun + 1;

    b = x + dx;
    b = x + 1;
    fb = trust_sub_secular_eqn(b, eigval, alpha, delta);
    itfun = itfun + 1;

    % find change of sign
    while ((fa > 0) == (fb > 0))

        dx = 2*dx;

        if ((fa > 0) ~= (fb > 0))
            break
        end
        b = x + dx;
        fb = trust_sub_secular_eqn(b, eigval, alpha, delta);
        itfun = itfun + 1;

        if (itfun > itbnd)
            break
        end
    end

    fc = fb;

    % main loop, exit from the middle of the loop
    while (fb ~= 0)
        % Ensure that b is the best result so far, a is the previous
        % value of b, and c is on hte opposit side of 0 from b
        if (fb > 0) == (fc > 0)
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        end

        if abs(fc) < abs(fb)
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        end

        % convergence test and possible exit
        if itfun > itbnd
            break
        end

        m = 0.5 * (c-b);
        rel_tol = 2.0 * tol * max(abs(b), 1.0);

        if (abs(m) <= rel_tol) || (abs(fb) < tol)
            break
        end

        % choose bisection or interpolation
        if (abs(e) < rel_tol) || (abs(fa) <= abs(fb))
            % bisection
            e = m;  d = e;
        else
            % interpolation
            s = fb/fa;
            if a == c
                % linear interpolation
                p = 2.0 * m * s;
                q = 1.0 - s;
            else
                % Inverse quad interpolation
                q = fa/fc;
                r = fb/fc;
                p = s * (2.0 * m * q * (q-r) - (b-a) * (r-1.0));
                q = (q-1.0) * (r-1.0) * (s-1.0);
            end
            
            if p > 0
                q = -q;
            else
                p = -p;
            end
            % Check if the interpolated point is acceptable
            if (2.0*p < 3.0*m*q - abs(rel_tol*q)) && (p < abs(0.5*e*q))
                e = d;
                d = p/q;
            else
                d = m;
                e = m;
            end
            %  end of iterpolation
        end

        % Next point
        a = b;
        fa = fb;
        if (abs(d) > rel_tol)
            b = b + d;
        else
            if b > c
                b = b - rel_tol;
            else
                b = b + rel_tol;
            end
        end

        fb = trust_sub_secular_eqn(b, eigval, alpha, delta);
        itfun = itfun + 1;
    end
end

