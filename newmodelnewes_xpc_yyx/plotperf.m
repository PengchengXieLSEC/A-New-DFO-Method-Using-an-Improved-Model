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
% clear
% clc

 load('T.mat');
 

  [np, ns] = size(T);
   fontsize = 30;
    linewidth = 2;
    penalty = 20;
    markersize = 10;

    logplot=1;

% Other colors, lines, and markers are easily possible:
colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];



% Compute ratios and divide by smallest element in each row.
r = T./repmat(min(T,[],2),1,ns);

% Replace all NaN's with twice the max_ratio and sort.
max_ratio = max(max(r));
r(isnan(r)) = 2*max_ratio;
r = sort(r);

% rtrans=r(:,5);
% r(:,5)=r(:,1);
% r(:,1)=rtrans;


% Plot stair graphs with markers.
hl = zeros(ns,1);
for s = 1:ns
    [xs,ys] = stairs(r(:,s),(1:np)/np);

    % Only plot one marker at the intercept
    if (xs(1)==1)
        vv = find(xs==1,1,'last');
        xs = xs(vv:end);   ys = ys(vv:end);
    end

    sl = mod(s-1,3) + 1; sc = mod(s-1,7) + 1; sm = mod(s-1,12) + 1;
    option1 = [char(lines(sl)) colors(sc) markers(sm)];
    if (logplot)
        hl(s) = semilogx(xs,ys,option1, 'LineWidth', linewidth);
    else
        hl(s) = plot(xs,ys,option1, 'LineWidth', linewidth);
    end
    hold on;
end

% Axis properties are set so that failures are not shown, but with the
% max_ratio data points shown. This highlights the "flatline" effect.
if (logplot) 
  axis([1 1.1*max_ratio 0 1]);
  twop = floor(log2(1.1*max_ratio));
  set(gca,'XTick',2.^[0:twop])
else
  % axis([1 1.1*max_ratio 0 1]);
end






hfig=figure(1);

% title('Performance profile');
 % titleString = sprintf('$\\tau = %.2e$', tauvalue);
% titleString= sprintf('$\\tau = 10^{-%d}$', 2*tauvalueindex-1);
 % title(titleString, 'Interpreter', 'latex');
 lgd = legend('Ours', 'Least Frob. norm', 'Powell', 'Least $H^2$ norm updating', 'Conn \& Toint', ...
     'NEWUOA', 'Fminsearch', 'Fminunc', 'CMA-ES', 'NMSMAX',...
     'Location', 'southeastoutside', 'Orientation', 'vertical', 'Interpreter', 'latex', 'fontsize', 20, 'Box', 'off', 'FontName', 'Times New Roman');

% lgd = legend('Ours', 'Least Frob. norm', 'Powell', 'Least $H^2$ norm updating', 'Conn \& Toint', ...
%     'NEWUOA', 'Fminsearch', 'Fminunc', 'CMA-ES', 'NMSMAX',...
%     'Location', 'southeastoutside', 'Orientation', 'vertical', 'Interpreter', 'latex', 'fontsize', fontsize, 'Box', 'off', 'FontName', 'Times New Roman');

% set(lgd, 'Orientation', 'horizontal', 'Location', 'southoutside'), 'Interpreter', 'latex', 'fontsize', 20, 'Box', 'off', 'FontName', 'Times New Roman';
% legend('boxoff'); % Turn off the legend box

% Arrange the legend in two rows
% numRows = 2;
% lgd.NumColumns = 5%ceil(lgd.NumEntries / numRows);




xlabel('$\alpha = \mathrm{NF}/\mathrm{NF}_{\min}$', 'fontsize', fontsize, 'interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\rho_a(\alpha)$', 'fontsize', fontsize, 'interpreter', 'latex', 'FontName', 'Times New Roman');
set(gca, 'FontSize', fontsize, 'LineWidth', linewidth, 'FontName', 'Times New Roman');
 % set(gcf, 'position', [100 100 700 1000]);
 % set(gcf, 'position', [100 100 1300 800]);
  % set(gcf, 'position', [100 100 1300 600]);
   set(gcf, 'position', [100 100 1400 600]);
  % set(gcf, 'position', [100 100 1700 1000]);
% set(hfig, 'visible', 'off');

 figname = strcat(int2str(int32(-log10(tauvalue))), '_perfwild');

 % figname = strcat(int2str(int32(-log10(tauvalue))), '_perfwild');
epsname = strcat(figname, '.eps');
saveas(hfig, epsname, 'epsc2');
resultsFolder = 'results';
% movefile(epsname, fullfile(resultsFolder, epsname));

% set(gcf, 'position', [100 100 1000 600]);
% Axis properties are set so that failures are not shown, but with the
% max_ratio data points shown. This highlights the "flatline" effect.
% axis([0 1.1*max_data 0 1]);







