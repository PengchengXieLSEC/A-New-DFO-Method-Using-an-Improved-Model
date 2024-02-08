function [H, g_hat, c_hat] = quad_frob_new(X_value, F_values, para, s)
    eps = 2.220446049250313e-16;
    tol_svd = eps.^5;
    % n = number of variables m = number of points
    [n, m] = size(X_value);

    H = zeros(n, n);
    g = zeros(n, 1);

    % Shift the points to the origin
    Y = X_value - diag(X_value(1:end, 1)) * ones(n, m);

    if (m < (n+1)*(n+2)/2)
        % 提取输入参数
        c_1 = para(1);
        c_2 = para(2);
        omega_1 = para(3);
        omega_2 = para(4);
        mu = para(5);
        theta = para(6);

        % 显式表达x0,xt,q
        x0 = X_value(:,1);
%         xk = x0;
        [~,xt_index] = max(F_values);
        xt = X_value(:,xt_index);
        q = xt - x0;

        % 计算K_G
        KG_tmp = c_2 * eye(n) - theta * (q * q');
        if cond(KG_tmp) < 1.0e+12
            KG = inv(KG_tmp);
        else
        % Make sure the condition number is not too high
            [U, S, VT] = svd(KG_tmp);
            indices = (S < tol_svd);
            S(indices) = tol_svd;
            Sinv = zeros(size(S,1), size(S,1));
            Sinv(~indices) = 1./S(~indices);
            % Get the coefficients
            KG = VT * Sinv * U';
        end

        index_nan = isnan(KG);
        KG(index_nan) = tol_svd;
        index_inf = isinf(KG);
        KG(index_inf) = 1e10;

        % 中间变量
        B = (2*c_1 + mu - theta) * eye(n) - theta^2 * norm(q,2)^2 * KG;
        A = zeros(m,m);
        X = zeros(n,m);%先算完再转置
        Yb = zeros(n,m);
        r = zeros(m,1);
        for i = 1:m
            X(:,i) = Y(:,i) + 0.5 * theta * KG' * (Y(:,i) * Y(:,i)') * q;
            Yb(:,i) = 0.5 * Y(:,i) + 0.25 * theta * KG * (Y(:,i) * Y(:,i)') * q;
            r(i) = F_values(i) - 0.5 * theta * omega_2 * Y(:,i)' * KG * (q * q') * Y(:,i);
            for j = 1:m
                A(i,j) = -0.125 * Y(:,i)' * KG * (Y(:,j) * Y(:,j)') * Y(:,i);
            end
        end
        y = (omega_2 * theta^2 * norm(q,2)^2 * KG * q) - (mu * omega_1 * s) + (theta * omega_2*q);
        
        % 右端向量
        b = [r; 0; y];
        % 系数矩阵
        W = [A, ones(m, 1), X'; ones(1, m), zeros(1,n+1); Yb, zeros(n,1), B];

        index_nan = isnan(W);
        W(index_nan) = tol_svd;
        index_inf = isinf(W);
        W(index_inf) = 1e10;
        
        [lambda_0] = quad_Frob_compute_coeffs(W, tol_svd, b, 'partial');

        % Grab the coeffs of linear terms (g) and the ones of quadratic terms
        % (H) for g.T s + s.T H s
        g = lambda_0(m+2:end);
        c = lambda_0(m+1);  % 新加入的c

        tmp = zeros(n,n);
        for i = 1:m
            tmp = tmp + lambda_0(i) * KG * (Y(:,i) * Y(:,i)');
        end
        H = theta * KG * g * q' - 0.25 * tmp + omega_2 * theta * KG * (q * q');
        
        g_hat = g - H * x0;
        c_hat = c + 0.5 * x0' * H * x0 - g' * x0;
        
    else  % Construct a full model
        % Here we have enough points. Solve the sys of equations.
        b = F_values;
        phi_Q = [];
        for i = 1:m
            y = Y(1:end, i);
            y = y(newaxis);    % turn y from 1D to a 2D array
            yguodu = (y.^2);
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
        g = lambda_0(1:n+1, 1:end);  % 1 n
        cont = n+1;
        H = zeros(n, n);

        for j = 1:n   
            H(j:n, j) = lambda_0(cont:cont + n - j, 1:end);  % n+1+n- 0  n-1 2n+1 n+2
            cont = cont + n - j;
        end

        H = H + H.T - diag(diag(H));
        
        
        c = [];
        g_hat = [];
        c_hat = [];
    end
end

