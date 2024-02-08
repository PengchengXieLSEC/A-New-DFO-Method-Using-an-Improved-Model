function [H, g_hat, c_hat] = quad_frob(X_value, F_values)
    eps = 2.220446049250313e-16;
    tol_svd = eps.^5;
    % n = number of variables m = number of points
    [n, m] = size(X_value);

    H = zeros(n, n);
    g = zeros(n, 1);

    % Shift the points to the origin
    Y = X_value - diag(X_value(1:end, 1)) * ones(n, m);

    if (m < (n+1)*(n+2)/2)
        
        b = [F_values; zeros(n+1, 1)];  % 图片的r
        A = 0.5 * (Y' * Y).^2;

       
        top = [A, ones(m, 1), Y'];
        temp = [ones(1, m),; Y];
        bottom = [temp, zeros(n+1, n+1)];
        W = [top; bottom];
        [lambda_0] = quad_Frob_compute_coeffs(W, tol_svd, b, 'partial');

        % Grab the coeffs of linear terms (g) and the ones of quadratic terms
        % (H) for g.T s + s.T H s
        g = lambda_0(m+2:end);
        
        c = lambda_0(m+1);  % 新加入的c

        H = zeros(n, n);
        for j = 1:m
            H = H + lambda_0(j) * reshape(Y(1:end, j),n, 1) * reshape(Y(1:end, j), 1, n);
        end
        
        g_hat = g - H * X_value(1:end, 1);
        c_hat = c + 0.5 * X_value(1:end, 1)' * H * X_value(1:end, 1) - g' * X_value(1:end, 1);
        
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

