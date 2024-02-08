% clear
% clc
% 
 load('T.mat');

tauvalueindex=2;
tauvalue=10^(-tauvalueindex);
  % T=Ttau(commonElements,:,tauvalueindex);

  [np, ns] = size(T);
   fontsize = 30;
    linewidth = 2;
    penalty = 20;
    markersize = 10;

% Other colors, lines, and markers are easily possible:
colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];
% Replace all NaN's with twice the max_ratio and sort.
max_data = max(max(T));
T(isnan(T)) = 2*max_data;
T = sort(T);

  filename = 'problems1';
    fileID = fopen(filename);
    C = textscan(fileID, '%s %f');
    tstn = C{2}(1:110);
    fclose(fileID);
     % tstn=tstn(commonElements);



% For each solver, plot stair graphs with markers.
hl = zeros(ns,1);
for s = 1:ns
    a=sort(T(:,s)./(tstn(:)+1));
    [xs,ys] = stairs(a,(1:np)/np);
    % [xs,ys] = stairs(T(:,s)./(tstn(mask)+1),(1:np)/np);
    sl = mod(s-1,3) + 1; sc = mod(s-1,7) + 1; sm = mod(s-1,12) + 1;
    option1 = [char(lines(sl)) colors(sc) markers(sm)];
    
    hl(s) = plot(xs,ys,option1,'LineWidth',linewidth);
    hold on;
end

  % axis([0 200 0 1]);

hfig=figure(1);

% title('Data profile');

titleString= sprintf('$\\tau = 10^{-%d}$', 2*tauvalueindex-1);
title(titleString, 'Interpreter', 'latex');
% 
legend('Ours', 'Least Frob. norm', 'Powell', 'Least $H^2$ norm updating', 'Conn \& Toint',...
    'NEWUOA', 'Fminsearch', 'Fminunc', 'CMA-ES', 'NMSMAX',...
    'Location', 'southeastoutside', 'Orientation', 'vertical', 'Interpreter', 'latex', 'fontsize', fontsize, 'Box', 'off', 'FontName', 'Times New Roman');
xlabel('$\beta = \mathrm{NF}/(n+1)$', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('$\delta_a(\beta)$', 'fontsize', fontsize, 'interpreter', 'latex');
set(gca, 'FontSize', fontsize, 'LineWidth', linewidth, 'FontName', 'Times New Roman');
% set(gcf, 'position', [100 100 1500 600]);
 % set(gcf, 'position', [100 100 700 1000]);
  % set(gcf, 'position', [100 100 700 1000]);
 set(gcf, 'position', [100 100 1400 600]);
% Axis properties are set so that failures are not shown, but with the
% max_ratio data points shown. This highlights the "flatline" effect.
  axis([0 100 0 1]);

% set(hfig, 'visible', 'off');

figname = strcat(int2str(int32(-log10(tauvalue))), '_datawild');
epsname = strcat(figname, '.eps');
saveas(hfig, epsname, 'epsc2');
resultsFolder = 'results';
% movefile(epsname, fullfile(resultsFolder, epsname));







