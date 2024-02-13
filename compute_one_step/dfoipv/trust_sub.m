function [s, val] = trust_sub(g, H, delta)
  %func_name可选rosen, sphere, ackley, beale, booth, bukin
  %   
    finfo_float_eps = 2.220446049250313e-16;
    tol = 10e-12;
    tol_seqeq = 10e-8;
    key = 0;
    itbnd = 50;
    lambda_0 = 0;
    s_factor = 0.8;
    b_factor = 1.2;
    n = size(g,1);
    coeff = zeros(n, 1);

    % convert H to full matrix if sparse
    T = sparse(H);
    H = full(T);

    % get the eigen value and vector
    [V, D] = eig(0.5 * (H' + H));
    count = 0;
    eigval = diag(D);
    % find the minimum eigen value
    [shunxu, index] = sort(eigval);
    jmin = index(1);
    mineig = shunxu(1);

    % depending on H, find a step size
    alpha = -V' * g;
    sig = sign(alpha(jmin)) + sum(alpha(jmin) == 0);
    sig = sig(1);

    % PSD case
    if mineig > 0
        lambda_0 = 0;
        coeff = alpha .* (1./eigval);
        s = V * coeff;
        % That is, s = -v (-v.T g./eigval)
        nrms = norm(s);
        if nrms < b_factor*delta
            key = 1;
        else
            laminit = [0];
        end
    else
        laminit = -mineig;
    end

    % Indefinite case
    if key == 0
        if trust_sub_secular_eqn(laminit, eigval, alpha, delta) > 0
          [b, c, count] = trust_sub_rfzero(laminit, itbnd, eigval, alpha, delta, tol);

          if abs(trust_sub_secular_eqn(b, eigval, alpha, delta)) <= tol_seqeq
              lambda_0 = b;
              key = 2;
              lam = lambda_0 * ones(n, 1);

              [coeff, s, nrms, w] = trust_sub_compute_step(alpha, eigval, coeff, V, lam);

              if (nrms > b_factor * delta || nrms < s_factor * delta)
                  key = 5;
                  lambda_0 = -mineig;
              end
          else
                key = 3;
                lambda_0 = -mineig;
          end
        else
            key = 4;
            lambda_0 = -mineig;
        end
        lam = lambda_0 * ones(n, 1);
        

        if key > 2
            arg = abs(eigval + lam) < 10 * (finfo_float_eps * max(abs(eigval), ones(n,1)));
            alpha(arg) = 0.0;
        end

        [coeff, s, nrms, w] = trust_sub_compute_step(alpha, eigval, coeff, V, lam);

        if key > 2 && nrms < s_factor * delta
            beta = sqrt(delta.^2 - nrms.^2);
            s = s + reshape(beta * sig * V(1:end, jmin), n, 1);
        end

        if key > 2 && nrms > b_factor * delta
            [b, c, count] = trust_sub_rfzero(laminit, itbnd, eigval, alpha, delta, tol);
            lambda_0 = b;
            lam = lambda_0 * ones(n, 1);

            [coeff, s, nrms, w] = trust_sub_compute_step(alpha, eigval, coeff, V, lam);
        end
    end

    % return the model prediction of the change in the objective with s
    val = g' * s + 0.5 * s' * H * s;
end

