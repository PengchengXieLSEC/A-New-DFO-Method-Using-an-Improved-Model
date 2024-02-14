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
function [H, g_hat, c_hat] = quad_frob_new(X_value, F_values, para, s, delta, eta, rho, HOLD)
  eps = 2.220446049250313e-16;
  tol_svd = eps .^ 5;
  % n = number of variables m = number of points
  [n, m] = size(X_value);

  H = zeros(n, n);
  g = zeros(n, 1);

  % Shift the points to the origin
  Y = X_value - diag(X_value(1:end, 1)) * ones(n, m);

  if (m < (n + 1) * (n + 2) / 2)

    % 显式表达x0,xt,q
    x0 = X_value(:, 1);
    %         xk = x0;
    [~, xt_index] = max(F_values);
    xt = X_value(:, xt_index);
    q =- (xt - x0);

 b1 = zeros(m, 1);  % 创建一个 m 行 1 列的零向量

% for i = 1:m
%     aterm = Y(:, i);
%     b1(i) = 0.5 * (aterm' * HOLD * aterm);
% end

% 预先计算 Y' 和 HOLD 的乘积
temp = Y' * HOLD; % temp 是一个 m x n 的矩阵

% 计算二次形式的向量化版本
b1 = 0.5 * sum(temp .* Y', 2);

    b = [F_values-b1; zeros(n + 1, 1)]; % 图片的r

% b = [F_values; zeros(n + 1, 1)]; % 图片的r


    A = 0.5 * (Y' * Y) .^ 2;

    top = [A, ones(m, 1), Y'];
    temp = [ones(1, m), ; Y];
    bottom = [temp, zeros(n + 1, n + 1)];
    W = [top; bottom];
    alpha = 1;
    beta = 1;
    if (norm(s, 2) < 1e-12) || (rho <= eta)
      alpha = 0;
      beta = 0;
    elseif delta == norm(s, 2)
      alpha = 0;
      beta = 1;
    elseif delta > norm(s, 2)
      alpha = 1;
      beta = 0;
    end
    %  alpha = 0;
    % beta = 0;
    p = s * s' / norm(s, 2) ^ 2;
    temp = -2 * alpha * eye(n) - 2 * beta * (eye(n) - p)' * (eye(n) - p);
    W(size(W, 1) - n + 1:end, size(W, 1) - n + 1:end) = temp;
    % [lambda_0] = quad_Frob_compute_coeffs(W, tol_svd, b, 'partial');
        lambda_0 = W \ b;
    % Grab the coeffs of linear terms (g) and the ones of quadratic terms
    % (H) for g.T s + s.T H s
    g = lambda_0(m + 2:end);

    c = lambda_0(m + 1); % 新加入的c

    H = zeros(n, n);
    % for j = 1:m
    %   H = H + lambda_0(j) * reshape(Y(1:end, j), n, 1) * reshape(Y(1:end, j), 1, n);
    % end


        Lambdaxpc =zeros(m,m);
Lambdaxpc = diag(lambda_0(1:m));

% 计算 H

H = Y * Lambdaxpc * Y';
H =H+HOLD;

    g_hat = g - H * X_value(1:end, 1);
    c_hat = c + 0.5 * X_value(1:end, 1)' * H * X_value(1:end, 1) - g' * X_value(1:end, 1);

  else % Construct a full model
    % Here we have enough points. Solve the sys of equations.
    b = F_values;
    phi_Q = [];
    for i = 1:m
      y = Y(1:end, i);
      y = y(newaxis); % turn y from 1D to a 2D array
      yguodu = (y .^ 2);
      aux_H = y * y.T - 0.5 * diag(yguodu(0));
      aux = [];
      for j = 1:n
        aux = [aux, aux_H(j:n, j)];
      end

      phi_Q = [phi_Q; aux];
    end

    W = [ones(m, 1), Y.T];
    W = [W, phi_Q];

    %lambda_0 = quad_Frob_compute_coeffs(W, tol_svd, b, option='full');

    % Retrieve the model coeffs (g) and (H)
    g = lambda_0(1:n + 1, 1:end); % 1 n
    cont = n + 1;
    H = zeros(n, n);

    for j = 1:n
      H(j:n, j) = lambda_0(cont:cont + n - j, 1:end); % n+1+n- 0  n-1 2n+1 n+2
      cont = cont + n - j;
    end

    H = H + H.T - diag(diag(H));

    c = [];
    g_hat = [];
    c_hat = [];
  end
end
