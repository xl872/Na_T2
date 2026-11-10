%% =====================================================================
%  Fig. 3D remake: ADC of the fast / slow diffusion ROI of 4 organoids, from
%  Fig3D_input_signals.csv (mean ROI signal per effective b value, 2 ROIs per organoid).
%
%  The fit reproduces the typed ADC values in Fig3D_input_sheetADC.csv exactly for Organoid 1 and 4; for Organoid 2 and 3
%         the typed values differ a little (printed below).
%
%  OUTPUTS (this folder)
%         Fig3D_ADC_table.csv / .xlsx       ADC (10^-4 mm^2/s, 95% CI) and R^2 of every fit (fast and slow of each organoid)
%         Fig3D_ADCfit_Organoid<n>.jpg/.eps  signal decay + fit, fast (blue) and slow (red) ROI, one figure per organoid
%         Fig3D_ADCfit_all.jpg/.eps         the 4 fits side by side
%         Fig3D_ADCbar.jpg/.eps             bar = mean, error bar = SEM, dots = organoids, paired t-test fast vs slow
%         Fig3D_ADCR2box.jpg/.eps           box plot of the fit R^2 (fast blue / slow red, as the bar plot), dots = organoids
%  The script changes to data/Fig3D itself.
%  =====================================================================
cd(fullfile(fileparts(mfilename('fullpath')),'..','data','Fig3D'));   % release layout: run in the dataset folder
clc; clear; close all;

%% ---------------- 0. settings ----------------------------------------------------
csvFile = 'Fig3D_input_signals.csv';  adcFile = 'Fig3D_input_sheetADC.csv';
% organoid name and ROI names; signals are the rows of csvFile for that organoid
D = struct('name',{'Organoid 1','Organoid 2','Organoid 3','Organoid 4'}, ...
           'roi',{{'ROI 1','ROI 2'},{'ROI 1','ROI 2'},{'left sample','right sample'},{'ROI 1','ROI 2'}});
excelADC = readmatrix(adcFile);   % typed ADC values (10^-4 mm^2/s), [fast slow] per organoid
useExcelADC = false;        % bar plot: false = the ADC fitted here (default, = the table), true = the values typed in the sheet
colF = [0 0 1];  colS = [1 0 0];     % fast = blue, slow = red (as Fig. 3D in Final Figures.pptx)
fs = 16;  axLW = 2;                  % font size (pt, Arial bold), axis line width
sigScale = 1e4;                      % signal axis in 10^4 a.u.
sigYLim = [0 2.5];  sigYTick = 0:0.5:2.5;   % same signal axis on all fit plots
bScale = 1e3;                        % b axis in 10^3 s/mm^2

%% ---------------- 1. fits ----------------------------------------------------------
ft  = fittype('A + C*exp(-x*ADC)', 'independent','x', 'coefficients',{'A','C','ADC'});
nO  = numel(D);  rows = {};  ADC = zeros(nO,2);  R2 = ADC;  F = cell(nO,2);  Bv = cell(nO,1);  Sv = cell(nO,2);
for o = 1:nO
    Tin = readtable(csvFile);  M = Tin{strcmp(Tin.organoid,D(o).name),{'b_s_mm2','roi_a','roi_b'}};  M = M(all(~isnan(M),2),:);
    b = M(:,1);  a = zeros(1,2);  r2 = a;  fo = cell(1,2);  ci = zeros(2,2);
    for r = 1:2
        y = M(:,r+1);
        [fo{r},g] = fit(b, y, ft, 'StartPoint',[0.8*min(y) max(y)-0.8*min(y) 5e-4], 'Lower',[0 0 0], 'Upper',[Inf Inf 0.01]);
        a(r) = fo{r}.ADC*1e4;  r2(r) = g.rsquare;  c = confint(fo{r});  ci(r,:) = c(:,3)'*1e4;      % 10^-4 mm^2/s
    end
    [~,iF] = max(a);  ord = [iF 3-iF];                                  % column 1 = fast, 2 = slow
    ADC(o,:) = a(ord);  R2(o,:) = r2(ord);  F(o,:) = fo(ord);  Bv{o} = b;  Sv(o,:) = {M(:,ord(1)+1), M(:,ord(2)+1)};
    cls = {'fast','slow'};
    for k = 1:2
        r = ord(k);
        rows(end+1,:) = {D(o).name, D(o).roi{r}, cls{k}, a(r), ci(r,1), ci(r,2), r2(r), excelADC(o,k), numel(b), min(b), max(b)}; %#ok<SAGROW>
    end
    fprintf('%s: fast = %s  ADC %.2f x 10^-4 mm^2/s (R^2 %.4f) | slow = %s  ADC %.2f (R^2 %.4f)   [sheet: %.2f / %.2f]\n', ...
            D(o).name, D(o).roi{ord(1)}, ADC(o,1), R2(o,1), D(o).roi{ord(2)}, ADC(o,2), R2(o,2), excelADC(o,:));
end
T = cell2table(rows, 'VariableNames',{'organoid','ROI','class','ADC_1e-4mm2s','ADC_CI95lo','ADC_CI95hi','R2','ADC_sheet_1e-4mm2s','n_b','b_min','b_max'});
writetable(T, 'Fig3D_ADC_table.csv');  writetable(T, 'Fig3D_ADC_table.xlsx', 'WriteMode','replacefile');
disp(T(:,[1:4 7 8]));
fprintf('mean +/- SD (n = %d):  fast ADC %.2f +/- %.2f, R^2 %.4f +/- %.4f | slow ADC %.2f +/- %.2f, R^2 %.4f +/- %.4f\n', ...
        nO, mean(ADC(:,1)), std(ADC(:,1)), mean(R2(:,1)), std(R2(:,1)), mean(ADC(:,2)), std(ADC(:,2)), mean(R2(:,2)), std(R2(:,2)));

%% ---------------- 2. fit figures (fast and slow ROI per organoid) --------------------
fAll = figure('Name','ADC fits','Position',[40 40 1600 400],'Color','w');  tl = tiledlayout(fAll,1,nO,'TileSpacing','compact','Padding','compact');
for o = 1:nO
    for pass = 1:2                                            % 1 = own figure, 2 = panel of the combined figure
        if pass == 1, fo1 = figure('Position',[100 100 430 360],'Color','w');  ax = axes(fo1);
        else,         ax = nexttile(tl,o); end
        hold(ax,'on');  bb = linspace(0, max(Bv{o})*1.05, 200)';  xs = bb/bScale;  h = gobjects(1,2);
        cols = [colF; colS];
        for k = 1:2
            plot(ax, xs, F{o,k}(bb)/sigScale, '-', 'Color',cols(k,:), 'LineWidth',2);
            h(k) = plot(ax, Bv{o}/bScale, Sv{o,k}/sigScale, 'o', 'MarkerFaceColor',cols(k,:), 'MarkerEdgeColor','none', 'MarkerSize',7);
        end
        xlim(ax,[0 xs(end)]);  ylim(ax,sigYLim);  yticks(ax,sigYTick);  xticks(ax,fewTicks([0 xs(end)],4));
        xlabel(ax,'b value (10^3 s/mm^2)');  ylabel(ax,'DWI signal (10^4 a.u.)');  title(ax,D(o).name);
        legend(h, {'Fast','Slow'}, 'Box','off', 'Location','northeast');
        styleAx(ax, axLW);
        if pass == 1
            fontsize(fo1,fs,"points");  boldArial(fo1);
            print(fo1, sprintf('Fig3D_ADCfit_Organoid%d.jpg',o), '-djpeg','-r300');  print(fo1, sprintf('Fig3D_ADCfit_Organoid%d.eps',o), '-depsc','-r300');
            close(fo1);
        end
    end
end
figure(fAll);  fontsize(fAll,fs,"points");  boldArial(fAll);
print(fAll, 'Fig3D_ADCfit_all.jpg', '-djpeg','-r300');  print(fAll, 'Fig3D_ADCfit_all.eps', '-depsc','-r300');

%% ---------------- 3. bar plot: fast vs slow ADC --------------------------------------
V = ADC;  if useExcelADC, V = excelADC; end
mu = mean(V);  se = std(V)/sqrt(nO);
[~,p,~,st] = ttest(V(:,1), V(:,2));                            % paired: fast and slow ROI of the same organoid
fprintf('bar plot (%s ADC): fast %.2f +/- %.2f SEM, slow %.2f +/- %.2f SEM; paired t-test t(%d) = %.2f, p = %.4f\n', ...
        ternary(useExcelADC,'sheet','fitted'), mu(1), se(1), mu(2), se(2), st.df, st.tstat, p);
fb = figure('Position',[100 100 340 340],'Color','w');  hold on   % same size as the R^2 box plot, so the text matches on the slide
bar(1, mu(1), 0.6, 'FaceColor','w', 'EdgeColor',colF, 'LineWidth',2);
bar(2, mu(2), 0.6, 'FaceColor','w', 'EdgeColor',colS, 'LineWidth',2);
errorbar(1:2, mu, se, 'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',10);
jit = 0.08*[-1.5 -0.5 0.5 1.5]';  jit = jit(1:nO);
scatter(1+jit, V(:,1), 60, colF, 'filled');  scatter(2+jit, V(:,2), 60, colS, 'filled');
yTop = ceil(max(V(:))*1.25);  yb = max(V(:))*1.1;
plot([1 1 2 2], [yb-0.15 yb yb yb-0.15], 'k-', 'LineWidth',1.5);
ht = text(1.5, yb, pstars(p), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
xlim([0.4 2.6]);  ylim([0 yTop]);  yticks(0:ternary(yTop>6,2,1):yTop);
xticks(1:2);  xticklabels({'Fast','Slow'});  ylabel('ADC (10^{-4} mm^2/s)');
fontsize(fb,fs,"points");  styleAx(gca, axLW);  boldArial(fb);
ht.FontSize = fs+6;                                           % significance stars
print(fb, 'Fig3D_ADCbar.jpg', '-djpeg','-r300');  print(fb, 'Fig3D_ADCbar.eps', '-depsc','-r300');

%% ---------------- 4. box plot of the fit R^2 (as R2bi_ROIbox_fig4 of 0518) ----------------------
fr = figure('Position',[100 100 340 340],'Color','w');
xg = [ones(nO,1); 2*ones(nO,1)];  yg = [R2(:,1); R2(:,2)];
hold on;  cols = [colF; colS];                               % fast = blue, slow = red, as the ADC bar plot
for k = 1:2
    boxchart(k*ones(nO,1), R2(:,k), 'BoxFaceColor','w', 'BoxFaceAlpha',1, 'BoxEdgeColor',cols(k,:), 'WhiskerLineColor',cols(k,:), 'MarkerStyle','none', 'BoxWidth',0.6, 'LineWidth',2);
    swarmchart(k*ones(nO,1), R2(:,k), 30, cols(k,:), 'filled', 'MarkerEdgeColor','none', 'XJitterWidth',0.35);
end
yl = [floor(min(yg)*100)/100 1];  ylim(yl);  yticks(fewTicks(yl,4));
xticks(1:2);  xticklabels({'Fast','Slow'});  xlim([0.4 2.6]);  ylabel('R^2');
fontsize(fr,fs,"points");  styleAx(gca, axLW);  boldArial(fr);
for k = 1:2, fprintf('%s: R^2 median = %.4f, 25-75%% = [%.4f %.4f], min %.4f\n', ternary(k==1,'fast','slow'), median(R2(:,k)), prctile(R2(:,k),[25 75]), min(R2(:,k))); end
print(fr, 'Fig3D_ADCR2box.jpg', '-djpeg','-r300');  print(fr, 'Fig3D_ADCR2box.eps', '-depsc','-r300');
fprintf('wrote Fig3D_ADC_table.csv/.xlsx, Fig3D_ADCfit_*.jpg/.eps, Fig3D_ADCbar.jpg/.eps, Fig3D_ADCR2box.jpg/.eps\n');

%% ================= local functions ===================================================
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
function s = pstars(p), if p < 0.001, s = '***'; elseif p < 0.01, s = '**'; elseif p < 0.05, s = '*'; else, s = 'n.s.'; end, end
function styleAx(ax, lw)                   % plot style of Final Figures.pptx: no grid, left + bottom axes only, thick black axes
    grid(ax,'off');  box(ax,'off');
    set(ax, 'XColor','k', 'YColor','k', 'LineWidth',lw, 'FontName','Arial', 'FontWeight','bold');
    set([ax.XLabel ax.YLabel ax.Title], 'Color','k', 'FontName','Arial', 'FontWeight','bold');
end
function boldArial(f)                      % all text of the figure black-bold Arial
    set(findall(f,'-property','FontName'),'FontName','Arial');  set(findall(f,'-property','FontWeight'),'FontWeight','bold');
    set(findall(f,'Type','legend'),'TextColor','k');
end
function t = fewTicks(lim, nMax)           % round tick values inside lim, at most nMax of them
    st = sort(reshape([1 2 2.5 5]'*10.^(-3:4), 1, []));
    for s = st, t = ceil(lim(1)/s-1e-9)*s : s : floor(lim(2)/s+1e-9)*s;  if numel(t) <= nMax, return; end, end
end
