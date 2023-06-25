%Codes for the paper entitled "A New Derivative-free Method Using an Improved Under-determined Quadratic Interpolation Model"
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

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

for i = 1:m
    aterm = Y(:, i);
    b1(i) = 0.5 * (aterm' * HOLD * aterm);
end

    b = [F_values-b1; zeros(n + 1, 1)]; % 图片的r
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
    p = s * s' / norm(s, 2) ^ 2;
    temp = -2 * alpha * eye(n) - 2 * beta * (eye(n) - p)' * (eye(n) - p);
    W(size(W, 1) - n + 1:end, size(W, 1) - n + 1:end) = temp;
    [lambda_0] = quad_Frob_compute_coeffs(W, tol_svd, b, 'partial');

    % Grab the coeffs of linear terms (g) and the ones of quadratic terms
    % (H) for g.T s + s.T H s
    g = lambda_0(m + 2:end);

    c = lambda_0(m + 1); % 新加入的c

    H = zeros(n, n);
    for j = 1:m
      H = H + lambda_0(j) * reshape(Y(1:end, j), n, 1) * reshape(Y(1:end, j), 1, n);
    end

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
