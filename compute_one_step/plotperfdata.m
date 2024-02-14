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
load('T.mat');
  [np, ns] = size(T);
  if (np == 0) 
    T = ones(1, ns);
    np = 1;
  end

  Tmin = min(T, [], 2);

    fontsize = 20;
    linewidth = 1;
    penalty = 20;
    markersize = 10;

    % Defining 10 different colors
colors = [0.01 0.01 0.01;         % Black
          0.99 0.01 0.01;         % Red
          0.01 0.99 0.01;         % Green
          0.01 0.01 0.99;         % Blue
          0.99 0.01 0.99;         % Magenta
          0.01 0.99 0.99;         % Cyan
          0.99 0.99 0.01;         % Yellow
          0.5 0.5 0.5;   % Gray
          0.99 0.5 0.01;       % Orange
          0.5 0.01 0.99];      % Purple
   % 
   % colors = [0.99 0 0;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.01 0.01 0.01;
   %            0.0 0.0 0.98];


markers = ['s', '^', '>', '+', 'x', 'o', '*', 'd', 'p', 'h'];
lines = {'-', '--', ':', '-.', '-', '--', ':', '-.', '--', ':'};
  

%%%%%%%%%%%%%%%%%Perf Profile%%%%%%%%%%%%%%%%%
r = zeros(np, ns);
for ip = 1:np
    r(ip, :) = T(ip, :) / Tmin(ip);
end
r = log2(r);
max_ratio = max(1.0e-6, max(max(r)));
r(isnan(r)) = penalty * max_ratio;
r = sort(r);

maker_indx = 1:20:2 * size(r, 1);

clf;
hfig = figure(1);
for is = 1:ns
    [xs, ys] = stairs(r(:, is), (1:np) / np);

    
    option = [char(lines(mod(is, 15) )), markers(mod(is, 15) ), colors(mod(is, 15) , :)];
    
    % optionpre = [char(linespre(mod(is, 15) + 1)), colors(mod(is, 15) + 1, :)];
    % indicespre = xs <= 0.02;
    % xspre = xs(indicespre);
    % yspre = ys(indicespre);
    % 
    % plot(xspre, yspre, optionpre, 'Linewidth', linewidth, 'Color', colors(mod(is, 15) + 1, :), 'HandleVisibility', 'off'); 
    % hold on;

    indices = xs ~= 0;
    xs_filtered = xs(indices);
    ys_filtered = ys(indices);
    plot(xs_filtered, ys_filtered, option, 'Linewidth', linewidth, 'MarkerIndices', maker_indx, 'MarkerSize', markersize, ...
        'Color', colors(mod(is, 15), :), 'DisplayName', 'Plot 2', 'HandleVisibility', 'on');
    hold on;
end

% axis([0 0.8 0.5 1]);
title('Performance profile', 'FontName', 'Times New Roman');
legend('Conn \& Toint', 'Least Frob. norm', 'Powell', 'Least $H^2$ norm updating', 'Ours',...
    'NEWUOA', 'Fminsearch', 'Fminunc', 'CMA-ES', 'NMSMAX',...
    'Location', 'southeast', 'Orientation', 'vertical', 'Interpreter', 'latex', 'fontsize', fontsize, 'Box', 'off', 'FontName', 'Times New Roman');
xlabel('$\log_2(\alpha), \quad \alpha = \mathrm{NF}/\mathrm{NF}_{\min}$', 'fontsize', fontsize, 'interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\rho_a(\alpha)$', 'fontsize', fontsize, 'interpreter', 'latex', 'FontName', 'Times New Roman');
set(gca, 'FontSize', fontsize, 'LineWidth', linewidth, 'FontName', 'Times New Roman');
set(gcf, 'position', [100 100 1000 600]);
set(hfig, 'visible', 'off');

figname = strcat(int2str(int32(-log10(tau))), '_perf');
epsname = strcat(figname, '.eps');
saveas(hfig, epsname, 'epsc2');
resultsFolder = 'results';
movefile(epsname, fullfile(resultsFolder, epsname));

hold off


%%%%%%%%%%%%%%%%%Data Profile%%%%%%%%%%%%%%%%%
% M = ceil(M / min(tstn(mask)));
M = ceil(M / min(tstn(:)));
D = NaN(ns, M);

for is = 1:ns
    for i = 1:M
        % D(is, i) = length(find(T(:, is) ./ tstn(mask) <= i)) / (mm - 1);
        D(is, i) = length(find(T(:, is) ./ tstn(:) <= i)) / (mm - 1);
    end
    for i = (M + 1):(1.2 * M)
        D(is, i) = D(is, M);
    end
end


hfig = figure(1);
for is = 1:ns
    if is == 1 || is == 2 || is == 3 || is == 4 || is == 5
        maker_indx = 1:1:20;
    end

    [xs, ys] = stairs([1:1.2 * M] / 1.2, D(is, :));
    option = [char(lines(mod(is, 15))), markers(mod(is, 15)), colors(mod(is, 15), :)];

    plot(xs, ys, option, 'Linewidth', linewidth, 'MarkerIndices', maker_indx, 'MarkerSize', markersize, 'Color', colors(mod(is, 15), :));
    hold on;
end

hold off;

% axis([0 15 0.45 0.75]);
title('Data profile');
legend('Conn \& Toint', 'Least Frob. norm', 'Powell', 'Least $H^2$ norm updating', 'Ours', 'Location', 'southeast', 'Orientation', 'vertical', 'Interpreter', 'latex', 'fontsize', fontsize, 'Box', 'off', 'FontName', 'Times New Roman');
xlabel('$\beta = \mathrm{NF}/(n+1)$', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('$\delta_a(\beta)$', 'fontsize', fontsize, 'interpreter', 'latex');
set(gca, 'FontSize', fontsize, 'LineWidth', linewidth, 'FontName', 'Times New Roman');
set(gcf, 'position', [100 100 1000 600]);
set(hfig, 'visible', 'off');

figname = strcat(int2str(int32(-log10(tau))), '_data');
epsname = strcat(figname, '.eps');
saveas(hfig, epsname, 'epsc2');
resultsFolder = 'results';
movefile(epsname, fullfile(resultsFolder, epsname));

