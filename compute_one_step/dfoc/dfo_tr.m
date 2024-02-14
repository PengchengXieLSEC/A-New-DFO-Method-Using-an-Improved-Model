function [res, iteration, f_hist] = dfo_tr(bb_func, x_initial, options, para)
  %UNTITLED7 此处显示有关此函数的摘要
  %   此处显示详细说明
  f_hist = [];

  % start timing and set the paramters
  t1 = clock;

  %%
  % see the param class for the description of parameters
all_options = struct("delta", 1.0, ... % initial delta (i.e. trust region radius)
    "tol_delta", 1e-10, ... % smallest delta to stop
    "max_delta", 100.0, ... % max possible delta
    "gamma1", 1, ... % radius shrink factor
    "max_iter", 1, ... % maximum number of iterations
    "eta0", 0.0, ... % step acceptance test (pred/ared) threshold
    "eta1", 0.25, ... 
    "eta2", 0.75, ... %eta1 in the paper
    "tol_f", 1e-15, ... % threshold for abs(fprev - fnew)- used to stop
    "gamma2", 1, ... % radius expansion factor
    'tol_norm_g', 1e-15, ... % threshold for norm of gradient- used to stop
    'tol_norm_H', 1e-10, ... % threshold for the (frobenius) norm of H
    "min_del_crit", 1e-8, ... % minimum delta for the criticality step
    "min_s_crit", 0.1); % minimum step for the criticality step};

  field = fieldnames(options);
  length = size(field, 1);
  for i = 1:length
    fieldname = string(field(i));
    if ismember(fieldname, fieldnames(all_options))
      all_options = setfield(all_options, fieldname, getfield(options, fieldname));
    end
  end

  delta = all_options.delta;
  tol_delta = all_options.tol_delta;
  max_delta = all_options.max_delta;
  gamma1 = all_options.gamma1;
  max_iter = all_options.max_iter;
  eta0 = all_options.eta0;
  eta1 = all_options.eta1;
  eta2 = all_options.eta2;
  tol_f = all_options.tol_f;
  gamma2 = all_options.gamma2;
  tol_norm_g = all_options.tol_norm_g;
  tol_norm_H = all_options.tol_norm_H;
  min_del_crit = all_options.min_del_crit;
  min_s_crit = all_options.min_s_crit;

  %%
  % set the verbosity parameter
  if ismember("verbosity", fieldnames(options))
    verb = options.verbosity;
    if verb == 0 || verb == 1 %or
      verbosity = verb;
    else
      error("verbosity option should be 0 or 1.")
    end
  else
    verbosity = 1;
  end
  % find the dimension of the problem
  % find the max and min number of points for the quadratic model
  n = size(x_initial, 1);
  % maxY = (n+1) * (n+2) / 2
maxY = n + 1;
  minY = n + 1;

  % iterations counters
  func_eval = 0; iter_suc = func_eval; iteration = iter_suc;

  % construct the intial sample set and evaluate the objective at them
  % by evaluating the blackbox function
  % x_initial = reshape(x_initial, n, 1);

  %%
  % func: _build_initial_sample,
  % input:x_initial, delta, options,
  % output:Y, nY
  func_n = size(x_initial, 1);
  %     if ~isempty(options) || ~ismember("sample_gen", fieldnames(options))
  if isempty(options) || ~ismember("sample_gen", fieldnames(options))
    option_build = "auto";
  else
    option_build = options.sample_gen;
  end

  if strcmp(option_build, "auto")
    Y = repmat(x_initial, 1, 2 * func_n) + 0.5 * delta * [eye(func_n), -eye(func_n)];
    Y = [x_initial, Y];
elseif strcmp(option_build, "manual")
    %         Y1 = reshape([150.00, 150.00, 230.00, 150.00], func_n, 1);
    %         Y2 = reshape([150.00, 8.00, 230.00, 250.00], func_n, 1);
    %         Y3 = reshape([150.00, 50.00, 230.00, 230.00], func_n, 1);
  Y1 = reshape([2.5 0], func_n, 1);
  Y2 = reshape([2 0.5], func_n, 1);
Y1=x_initial+reshape([1 0], func_n, 1);
Y2=x_initial+reshape([0 1], func_n, 1);
  %Y3 = reshape([2.5 0], func_n, 1);
  Y = [x_initial, Y1, Y2];
  else
    error("The sample option should be manual or auto.")
  end

  nY = size(Y, 2);

  %%
  f_values = zeros(nY, 1);

  for i = 1:nY
    f_values(i) = bb_func(reshape(Y(1:end, i), n, 1));
    func_eval = func_eval + 1;
  end

  % find the point with minimum objective. set it as the first center
  [~, ind_f_sort] = min(f_values(1:end, 1));
  x = reshape(Y(:, ind_f_sort), n, 1);
  f = f_values(ind_f_sort, 1);

  % print out the initial evaluation results
  %if verbosity
    fprintf ("\n Iteration Report \n")
    fprintf ('|it |suc|   obj  | TR_radius |  rho  | |Y|  \n')
    fprintf ("| %d |---| %0.4f |   %0.3f   | ----- | %d \n", iteration, f, delta, nY)
  %end

  % Start the TR main loop
  H_hist = zeros(n, n, 1);
  g_hist = zeros(n, 1, 1);
  c_hist = [0];
  g_hat_hist = zeros(n, 1, 1);
  c_hat_hist = [0];
  % Q_value = zeros(nY, 1);

  s = [0; 0]; % s = xk - xk_1
  step = 1;
  while 1
    success = 0;

    % shift_Y = Y - np.dot(np.diag(Y[:, 0]), np.ones((n, nY)))
    % print('被用来减去的Q_value的Y', Y)
    % for j = 1:nY
      % Yzhuan = reshape(Y(1:end, j), 1, n);
      % g_hat_hist_rot = reshape(g_hat_hist(1:end, 1, end), 1, n);
      % Q_value(j) = 0.5 * Yzhuan * H_hist(1:end, 1:end, end) * Y(1:end, j) + g_hat_hist_rot * Y(1:end, j) + c_hat_hist(end);
    % end
    % print('被用来减去的Q_value', Q_value)
    % print('被用来减去的f_value', f_values)

    % f_min_Q_values = f_values - Q_value;
    % print('f_min_Q', f_min_Q_values)
    % construct the quadratic model
    if step == 1
      [H, g_hat, c_hat] = quad_frob(Y, f_values);
    else
      %[H, g_hat, c_hat] = quad_frob(Y, f_values);
      [H, g_hat, c_hat] = quad_frob_new(Y, f_values, para, s, delta, eta2, rho);
    end
    % 加前一步的H和g，没有c; 第k轮，加上第k-1的H和g，c  *** 2
    % print('Hg',(H, g))
    % Q_test = zeros(nY, 1);
    % for j = 1:nY
    %   Q_test(j) = 0.5 * Y(1:end, j)' * H * Y(1:end, j) + g_hat' * Y(1:end, j) + c_hat;
    % end

    %H = H + H_hist(1:end, 1:end, end); H_hist(1:end, 1:end, end+1)=H;

    %g_hat = g_hat + g_hat_hist(1:end, 1, end); g_hat_hist(1:end, 1, end+1)=g_hat;
    %c_hat = c_hat + c_hat_hist(end); c_hat_hist(end+1)=c_hat;
    guodu = H * Y(1:end, 1);
    g = g_hat + guodu(:);
    c = c_hat - 0.5 * Y(1:end, 1)' * H * Y(1:end, 1) - g' * Y(1:end, 1);

    normg = norm(g);

    if normg <= tol_norm_g || delta < tol_delta
      break
    end

    if norm(H, 'fro') > tol_norm_H % make sure the trust can work
      [s, val] = trust_sub(g, H, delta); % TARGET DELTA
    else % otherwise take a steepest descent step less than delta
      s =- (delta / normg) * g;
      val = g' * s + 0.5 * s' * H * s;
    end

    fmod = f + val;

    if abs(val) < tol_f || iteration > max_iter
      break
    end

    xtrial = x + s;

    ftrue = bb_func(xtrial);
    func_eval = func_eval + 1;

    pred = f - fmod;
    ared = f - ftrue;
    rho = ared / pred;

    if rho >= eta0 && f - ftrue > 0

      success = 1;
      iter_suc = iter_suc + 1;

      x = xtrial;
      f = ftrue;
      if rho >= eta2 % If confident, increase delta
        delta = min(gamma2 * delta, max_delta);
      end
    else
      if nY >= minY
        % decrease the TR radius
        delta = gamma1 * delta;
      end
    end

    iteration = iteration + 1;

    if size(f_hist, 1) == 0
      f_hist = f;
    else
      f_hist(end + 1, end) = f;
    end

    %if verbosity
      fprintf ("| %d | %d | %.4f |   %0.3f   | %0.3f | %d \n", iteration, success, f, delta, rho(1, 1), nY)
      fprintf ("%.4f \n", x');
    %end

    Ynorms = Y - repmat(x, 1, nY);
    Ynorms = sqrt(sum(Ynorms .^ 2, 1));
    [Ynorms, index] = sort(Ynorms);

    Y = Y(1:end, index);
    f_values = f_values(index);

    if success
      reached_max = 1;
      if nY < maxY
        reached_max = 0;
      end

      if reached_max % substitute
        Y(1:end, nY) = xtrial';
        f_values(nY) = ftrue;
      else % add
        nY = nY + 1;
        Y = [Y, xtrial];
        f_values = [f_values; ftrue];
      end
    else % if not successful
      if nY >= maxY
        if norm(xtrial - x) <= Ynorms(nY)
          Y(1:end, nY) = xtrial';
          f_values(nY) = ftrue;
        end
      else
        nY = nY + 1;
        Y = [Y, xtrial];
        f_values = [f_values; ftrue];
      end
    end

    Ynorms = Y - repmat(x, 1, nY);
    Ynorms = sqrt(sum(Ynorms .^ 2, 1));
    [Ynorms, index] = sort(Ynorms);

    Y = Y(1:end, index);
    f_values = f_values(index);

    if delta < min_del_crit && norm(s) < min_s_crit
      selection = (Ynorms < delta * 100);
      Y = Y(1:end, find(selection, 1, "first"));
      f_values = f_values(selection);
      nY = size(Y, 2);
    end

    % 从第二次开始，每次执行新quad_frob_new
    step = 2;
  end

  t2 = clock;
  time_consump = etime(t2, t1);
  if verbosity
    fprintf ('*****************REPORT************************\n');
    fprintf ("Total time is %.2f seconds.\n", time_consump);
    fprintf ("Norm of the gradient of the model is %.2f.\n", normg);
    fprintf ('***************Final Report**************\n');
    fprintf ("| iter | success | fevals | final fvalue | final tr_radius|\n");
    fprintf ("|  %d  |    %d   |  %.2f  |      %.2f    |  %.2f  \n", iteration, iter_suc, func_eval, f, delta)
  end
  res = struct("x", x, "fun", f, "iteration", iteration, "iter_suc", iter_suc, "func_eval", func_eval, "delta", delta);
end
