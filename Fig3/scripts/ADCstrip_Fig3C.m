%% =====================================================================
%  Fig. 3C panel (Fig. 3C of the paper):  per organoid  [ADC map] [ADC on the anatomy + strip ROI] [signal decay + ADC fit per cell]
%
%  The masks are the saved datasets ROI_ADCstrip+orig in each dataset folder (1 / 2 / 3 = cell 1 / 2 / 3).
%
%  The script changes to the release data folder itself.
%
%  ROWS   Organoid #1 : DWI_scan613 (1H scan 613)
%         Organoid #2 : org2set = '627' -> DWI_scan627 (scan 627, default)
%                                    org2set = '175' -> DWI_scan175    (scan 175)
%  SLICES are given the way AFNI shows them: org(o).view = 'sagittal' | 'coronal' | 'axial' and org(o).slice = the slice
%         number of that AFNI image viewer (0-based).  The dataset orientation (RSA here: i = R->L, j = S->I, k = A->P)
%         is read from the header and the panel is drawn like the AFNI viewer draws it:
%           sagittal : across = A->P, down = S->I, slices counted L->R  -> 613New sagittal 31 = array index i = 80-31 = 49
%           coronal  : across = R->L, down = S->I, slices counted A->P  -> 627    coronal  40 = array index k = 41
%
%  OUTPUTS (200_proc)      ADCstrip_Fig3C<tag>.jpg / .eps, ADCstrip_Fig3C<tag>_cells.csv (fit per cell),
%                          ADCstrip_Fig3C<tag>_decay.csv (mean, SE, n per cell and b value)      (<tag> = '' for 627 | '_175')
%  =====================================================================
cd(fullfile(fileparts(mfilename('fullpath')),'..','data'));   % release layout: run in the dataset folder
clc; clear; close all;

%% ---------------- 0. settings ----------------------------------------------------
base = pwd;                                          % .../data
RAW  = fullfile('bruker','20250801_201306_Mouse_23Naplus1_08012025_1_81');   % release layout
org2set = '627';            % dataset of Organoid #2: '627' = DWI_scan627 (default), '175' = DWI_scan175

org(1).label = 'Organoid #1';  org(1).folder = 'DWI_scan613';  org(1).scan = 613;
org(1).view  = 'sagittal';     org(1).slice  = 31;                    % AFNI sagittal viewer, slice 31
org(1).clim  = [];                                                    % ADC colour range (um^2/s); [] = 2nd-98th percentile of the ROI core
switch org2set
    case '627', org(2).folder = 'DWI_scan627';  org(2).scan = 627;  org(2).view = 'coronal';  org(2).slice = 40;  outTag = '';
    case '175', org(2).folder = 'DWI_scan175';     org(2).scan = 175;  org(2).view = 'coronal';  org(2).slice = 19;  outTag = '_175';
    otherwise,  error('org2set must be ''175'' or ''627''');
end
org(2).label = 'Organoid #2';  org(2).clim = [];

cellMM    = 0.3;
nCells    = 3;
mapModel  = 'fit';          % 'fit'  : ADC map of the displayed slice fitted here, voxel by voxel
                            % 'file' : ADC+orig of myADCProt.m
bMax      = Inf;            % use only b values <= bMax in the fits
smoothSigma = 0.5;          % 3-D Gaussian smoothing of every b volume before the fits, as in myADCProt.m (0 = none).
                            % Neighbouring voxels share signal after smoothing, so the SE bars are somewhat too small.
axisDir   = 'ij';           % 'ij' = panel drawn as the AFNI viewer draws it (published Fig. 3C); 'xy' = upside down
marginMM  = 0.5;            % field of view: organoid ROI box of the slice stack + this margin (mm)
ovAlpha   = 0.6;            % opacity of the ADC colours on the anatomy (middle column)
showCellNumbers = false;    % true: write 1 2 3 beside the cells on the map (the dashed outline colours = the curve colours either way)
col = [0.165 0.471 0.839;   % blue    fastest cell  (colorByADC = true)  | cell 1 (colorByADC = false)
       0.106 0.686 0.478;   % green   middle                             | cell 2
       0.922 0.408 0.204];  % orange  slowest                            | cell 3
colorByADC = true;          % true: colours follow the fitted ADC (fast -> slow = blue - green - orange) on the map and the
                            % decay plot; false: colours follow the position along the strip (cell 1 / 2 / 3)
mk  = {'o','s','^'};        % marker per cell (identity never by colour alone)
ink = [0 0 0];              % text / axis colour (black, as in Final Figures.pptx)
fs  = 20;                   % font size (pt) of all text; Arial bold as in Final Figures.pptx
axLW = 2;                   % axis line width
adcScale = 100;             % ADC shown in 10^-4 mm^2/s (= um^2/s / 100)
adcLab = 'ADC (10^{-4} mm^2/s)';
sigScale = 1e4;             % decay plot: signal in 10^4 a.u.
sigYLim = [0 4];  sigYTick = 0:1:4;   % y axis of the decay plots (10^4 a.u.)
bScale = 1e3;               % decay plot: b in 10^3 s/mm^2

%% ---------------- 1. load, ADC map of the slice, place the strip, fit the cells --------
ftS = fittype('A + C*exp(-x*ADC)', 'independent','x', 'coefficients',{'A','C','ADC'});
fitS = @(b,y) fit(b, y, ftS, 'StartPoint',[0.8*min(y) max(y)-0.8*min(y) 5e-4], 'Lower',[0 0 0], 'Upper',[Inf Inf 0.01]);
rows = {};  drow = {};
for o = 1:numel(org)
    d = fullfile(base, org(o).folder);
    [adc,Info] = BrikLoad(fullfile(d,'ADC+orig.BRIK'));  adc = double(adc);  ND = size(adc);   % myADCProt.m map, 1e-6 mm^2/s = um^2/s
    roi  = BrikLoad(fullfile(d,'ROI+orig.BRIK')) > 0;
    anat = double(BrikLoad(fullfile(d,'allave+orig.BRIK')));
    dat  = double(BrikLoad(fullfile(d,'allave_t+orig.BRIK')));                                  % b-value series
    res  = abs(Info.DELTA(1));  B = readBvalues(fullfile(RAW,num2str(org(o).scan),'method'));
    assert(numel(B) == size(dat,4), '%s: %d b values in the method file, %d volumes in allave_t', org(o).folder, numel(B), size(dat,4));
    if smoothSigma > 0, for t = 1:size(dat,4), dat(:,:,:,t) = imgaussfilt3(dat(:,:,:,t), smoothSigma); end, end

    % ---- into the layout of the AFNI viewer: V(x across, y down, slice number + 1)
    [perm, flp] = afniView(Info.ORIENT_SPECIFIC, org(o).view);
    adc = toView(adc,perm,flp);  roi = toView(roi,perm,flp);  anat = toView(anat,perm,flp);  dat = toView(dat,perm,flp);
    kk = org(o).slice + 1;  w = max(1, round(cellMM/res));  use = B <= bMax;
    assert(kk >= 1 && kk <= size(roi,3) && any(any(roi(:,:,kk))), '%s: %s slice %d does not cut the organoid ROI', org(o).label, org(o).view, org(o).slice);
    arrIdx = ternary(flp(3), size(roi,3) - org(o).slice, kk);
    fprintf('\n%s = %s (scan %d): %dx%dx%d @ %.2f mm, AFNI %s slice %d = array index %d along dataset axis %d (ROI %d voxels there), cell %dx%d voxels, strip %.2f x %.2f mm\n', ...
            org(o).label, org(o).folder, org(o).scan, ND, res, org(o).view, org(o).slice, arrIdx, perm(3), nnz(roi(:,:,kk)), w, w, w*res, nCells*w*res);

    % ---- ADC map of the slice
    M2 = roi(:,:,kk);
    if strcmpi(mapModel,'fit')
        A2 = zeros(size(M2));  R2 = A2;  [vx,vy] = find(M2);
        for q = 1:numel(vx)
            y = squeeze(dat(vx(q),vy(q),kk,:));  [f1,g1] = fitS(B(use), y(use));
            A2(vx(q),vy(q)) = f1.ADC*1e6;  R2(vx(q),vy(q)) = g1.rsquare;
        end
        fprintf('  ADC map: A + C*exp(-b*ADC) in %d voxels, median %.0f um^2/s, R^2 median %.3f (min %.3f), %d voxel(s) on the ADC bound\n', ...
                numel(vx), median(A2(M2)), median(R2(M2)), min(R2(M2)), nnz(A2(M2) >= 0.01*1e6*0.999 | A2(M2) <= 0));
    else
        A2 = adc(:,:,kk);  R2 = nan(size(A2));  fprintf('  ADC map: ADC+orig of myADCProt.m (C*exp(-b*ADC), no plateau)\n');
    end
    lim = org(o).clim;
    if isempty(lim)
        core = imerode(M2, strel('disk',1,0));  if nnz(core) < 10, core = M2; end
        lim = [floor(prctile(A2(core),2)/50) ceil(prctile(A2(core),98)/50)]*50;  fprintf('  ADC colour range [%g %g] um^2/s\n', lim);
    end

    % ---- strip: label image L (0 = outside, 1..nCells = cell number) from the saved mask ROI_ADCstrip+orig
    Ls = toView(BrikLoad(fullfile(d,'ROI_ADCstrip+orig.BRIK')),perm,flp);  L = Ls(:,:,kk);

    % ---- signal decay (mean, SE over the voxels of the cell) and fit per cell
    S = zeros(nCells, numel(B));  SE = S;  cellADCmap = zeros(nCells,1);  cellADCsd = cellADCmap;  cellR2 = cellADCmap;  cellR2sd = cellADCmap;  nv = zeros(nCells,1);  fitA = zeros(nCells,1);  fo = cell(nCells,1);
    for c = 1:nCells
        m = L == c;  nv(c) = nnz(m);
        for t = 1:numel(B), v = dat(:,:,kk,t);  S(c,t) = mean(v(m));  SE(c,t) = std(v(m))/sqrt(nv(c)); end
        cellADCmap(c) = mean(A2(m));  cellADCsd(c) = std(A2(m));              % voxel-wise ADCs of the map inside the cell
        cellR2(c) = mean(R2(m));  cellR2sd(c) = std(R2(m));                  % R^2 of those voxel-wise fits (NaN with mapModel = 'file')
        [fo{c},g] = fitS(B(use), S(c,use)');  fitA(c) = fo{c}.ADC*1e6;  ci = confint(fo{c})*1e6;
        fprintf('  cell %d: %d voxels, ADC map mean %6.1f | cell fit ADC = %6.1f um^2/s = %.2f x 10^-4 mm^2/s (95%% CI %.0f-%.0f), A = %.0f, C = %.0f, plateau A/(A+C) = %.0f%%, R^2 %.4f\n', ...
                c, nv(c), cellADCmap(c), fitA(c), fitA(c)/adcScale, ci(1,3), ci(2,3), fo{c}.A, fo{c}.C, 100*fo{c}.A/(fo{c}.A+fo{c}.C), g.rsquare);
        rows(end+1,:) = {org(o).label, org(o).folder, org(o).scan, org(o).view, org(o).slice, c, nv(c), cellADCmap(c), cellADCsd(c), cellR2(c), cellR2sd(c), fitA(c), fitA(c)/adcScale, ci(1,3), ci(2,3), fo{c}.A, fo{c}.C, g.rsquare}; %#ok<SAGROW>
        for t = 1:numel(B), drow(end+1,:) = {org(o).label, c, B(t), S(c,t), SE(c,t), nv(c), use(t)}; end %#ok<SAGROW>
    end
    fprintf('  ADC along the strip: %s um^2/s (%s)\n', sprintf('%.0f  ', fitA), ternary(all(diff(fitA) < 0), 'falls from cell 1 to 3', 'not monotonic'));
    colIdx = 1:nCells;                                                  % colour of each cell: row of col
    if colorByADC, [~,rk] = sort(fitA,'descend');  colIdx(rk) = 1:nCells; end
    fprintf('  colours (cell 1..%d): %s\n', nCells, strjoin(arrayfun(@(c) sprintf('cell %d = %s', c, ternary(colIdx(c)==1,'blue',ternary(colIdx(c)==2,'green','orange'))), 1:nCells, 'UniformOutput',false), ', '));
    fprintf('  ADC of the voxels in each cell (%s), um^2/s:\n', ternary(strcmpi(mapModel,'fit'), 'voxel-wise A + C*exp(-b*ADC) fits', 'ADC+orig of myADCProt.m'));
    for c = 1:nCells, fprintf('    cell %d:  mean %6.1f   std %6.1f   (n = %d voxels, min %.0f, max %.0f)\n', c, cellADCmap(c), cellADCsd(c), nv(c), min(A2(L == c)), max(A2(L == c))); end
    fprintf('  R^2 of the voxel-wise fits in each cell:\n');
    for c = 1:nCells, fprintf('    cell %d:  mean %.3f   std %.3f   (n = %d voxels, min %.3f, max %.3f)\n', c, cellR2(c), cellR2sd(c), nv(c), min(R2(L == c)), max(R2(L == c))); end

    R(o) = struct('adc',A2, 'mask',M2, 'anat',anat(:,:,kk), 'anatHi',prctile(reshape(anat(:,:,kk),[],1),99.5), 'L',L, 'res',res, 'B',B, 'use',use, ...
                  'S',S, 'SE',SE, 'fo',{fo}, 'fitA',fitA, 'colIdx',colIdx, 'lim',lim, 'box',roiBox(roi, round(marginMM/res))); %#ok<SAGROW>
end
T = cell2table(rows, 'VariableNames', {'organoid','folder','scan','afniView','afniSlice','cell','nVoxels','ADCmap_mean_um2s','ADCmap_std_um2s','ADCmap_R2_mean','ADCmap_R2_std', ...
               'ADCfit_um2s','ADCfit_1e-4mm2s','ADCfit_CI95lo','ADCfit_CI95hi','A_plateau','C_amplitude','R2'});
writetable(T, fullfile(base,['ADCstrip_Fig3C' outTag '_cells.csv']));
writetable(cell2table(drow, 'VariableNames', {'organoid','cell','b_s_mm2','meanSignal','SE','nVoxels','usedInFit'}), fullfile(base,['ADCstrip_Fig3C' outTag '_decay.csv']));
fprintf('\n  wrote ADCstrip_Fig3C%s_cells.csv and ADCstrip_Fig3C%s_decay.csv\n', outTag, outTag);

%% ---------------- 2. figure ------------------------------------------------------------
nO = numel(org);  J = jet(256);
figure('Name','ADC strip (Fig. 3C style)','Position',[40 40 1250 400*nO],'Color','w');
for o = 1:nO
    r = R(o);  lim = r.lim;
    cidx = min(256, max(1, round((r.adc - lim(1))/diff(lim)*255) + 1));  C = reshape(J(cidx(:),:), [size(r.adc) 3]);
    G = repmat(min(1, r.anat/r.anatHi), 1, 1, 3);  M3 = repmat(r.mask, 1, 1, 3);
    P1 = zeros(size(C));  P1(M3) = C(M3);                                 % ADC map on black
    P2 = G;  P2(M3) = (1-ovAlpha)*G(M3) + ovAlpha*C(M3);                  % ADC colours blended into the anatomy (opaque -> prints as EPS)

    % -- ADC map
    y0 = 1 - o/nO + 0.22/nO;  h0 = 0.66/nO;                               % row geometry (normalised figure units)
    ax = axes('Position',[0.09 y0 0.24 h0]);  image(permute(P1,[2 1 3]));  axis image off;  set(ax,'YDir',ydir(axisDir));
    xlim(r.box(1:2));  ylim(r.box(3:4));  colormap(ax,J);  clim(ax,lim);
    % hold on;  scaleBar(r.box, 1/r.res, axisDir);                       % 1 mm bar: added in the pptx instead
    hb = colorbar(ax,'Location','manual','Position',[0.065 y0+0.1*h0 0.011 0.8*h0],'AxisLocation','out');  hb.YAxisLocation = 'left';
    hb.Ticks = ceil(lim(1)/adcScale)*adcScale : adcScale : floor(lim(2)/adcScale)*adcScale;   % whole 10^-4 mm^2/s only
    hb.TickLabels = compose('%g', hb.Ticks/adcScale);                                    % um^2/s -> 10^-4 mm^2/s
    ylabel(hb, adcLab, 'FontWeight','bold');  hb.Color = ink;  hb.LineWidth = axLW;  hb.FontName = 'Arial';  hb.FontWeight = 'bold';
    % title(sprintf('%s: ADC, %s slice %d', org(o).label, org(o).view, org(o).slice), 'FontWeight','normal', 'Color',ink);

    % -- ADC on the anatomy + strip + scale bar
    ax = axes('Position',[0.37 y0 0.24 h0]);  image(permute(P2,[2 1 3]));  axis image off;  set(ax,'YDir',ydir(axisDir));  hold on
    xlim(r.box(1:2));  ylim(r.box(3:4));
    for c = 1:nCells
        [ci,cj] = find(r.L == c);  pos = [min(ci)-0.5 min(cj)-0.5 max(ci)-min(ci)+1 max(cj)-min(cj)+1];
        rectangle('Position',pos, 'EdgeColor','w', 'LineWidth',2.6);                       % white underlay: visible on any colour
        rectangle('Position',pos, 'EdgeColor',col(r.colIdx(c),:), 'LineWidth',1.6, 'LineStyle','--');
    end
    if showCellNumbers, labelCells(r.L, nCells, r.box); end
    % scaleBar(r.box, 1/r.res, axisDir);   (the 1 mm bar is on the ADC map, left)
    % title('ADC on anatomy, strip ROI', 'FontWeight','normal', 'Color',ink);

    % -- decay (mean +/- SE over the voxels of the cell) + fit per cell
    ax = axes('Position',[0.70 y0 0.27 h0]);  hold on;  bb = linspace(0, max(r.B(r.use))*1.2, 200)';  hl = gobjects(nCells,1);
    for c = 1:nCells
        cc = col(r.colIdx(c),:);
        hl(c) =plot(bb/bScale, r.fo{c}(bb)/sigScale, '-', 'Color',cc, 'LineWidth',2);
        errorbar(r.B/bScale, r.S(c,:)/sigScale, r.SE(c,:)/sigScale, 'LineStyle','none', 'Color',cc, 'LineWidth',2, 'CapSize',8);
        % hl(c) = plot(r.B(r.use)/bScale, r.S(c,r.use)/sigScale, mk{c}, 'MarkerSize',2, 'MarkerFaceColor',cc, 'MarkerEdgeColor','w', 'LineWidth',1);
        if any(~r.use), plot(r.B(~r.use)/bScale, r.S(c,~r.use)/sigScale, mk{c}, 'MarkerSize',8, 'MarkerFaceColor','w', 'MarkerEdgeColor',cc, 'LineWidth',1); end   % not fitted: open
    end
    xlim([0 bb(end)/bScale]);  ylim(sigYLim);  yticks(sigYTick);  set(ax,'TickDir','in','Layer','top');
    xlabel('b value (10^3 s/mm^2)');  ylabel('DWI signal (10^4 a.u.)');     % no legend: the ADC values are in the table
    styleAx(ax, ink, axLW);
end
fontsize(gcf,fs,"points")
set(findall(gcf,'-property','FontName'),'FontName','Arial');  set(findall(gcf,'-property','FontWeight'),'FontWeight','bold');
for a = findall(gcf,'Type','axes')', a.Toolbar.Visible = 'off'; end      % otherwise the axes toolbar icons can end up in the exported image
print(gcf, fullfile(base,['ADCstrip_Fig3C' outTag '.jpg']), '-djpeg', '-r300');
print(gcf, fullfile(base,['ADCstrip_Fig3C' outTag '.eps']), '-depsc', '-r300');
fprintf('  wrote ADCstrip_Fig3C%s.jpg / .eps\n', outTag);


%% ================= local functions ===================================================
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
function styleAx(ax, ink, lw)              % plot style of Final Figures.pptx: no grid, left + bottom axes only, thick black axes, black bold Arial
    grid(ax,'off');  box(ax,'off');
    set(ax, 'XColor',ink, 'YColor',ink, 'LineWidth',lw, 'FontName','Arial', 'FontWeight','bold');
    set([ax.XLabel ax.YLabel ax.Title], 'Color',ink, 'FontName','Arial', 'FontWeight','bold');
end
function d = ydir(axisDir), if strcmpi(axisDir,'xy'), d = 'normal'; else, d = 'reverse'; end, end
function B = readBvalues(methodFile)       % effective b values (s/mm^2) from the Bruker method file, as in myADCProt.m
    lines = regexp(fileread(methodFile), '\r?\n', 'split');
    i1 = find(strncmp(lines,'##$PVM_DwEffBval',16), 1);  assert(~isempty(i1), 'no PVM_DwEffBval in %s', methodFile);
    i2 = i1 + find(strncmp(lines(i1+1:end),'##',2) | strncmp(lines(i1+1:end),'$$',2), 1);
    B = sscanf(strjoin(lines(i1+1:i2-1),' '), '%f');
end
function [perm, flp] = afniView(orientCodes, view)   % dataset axes -> [across, down, slice order] of an AFNI image viewer
    % AFNI orientation codes (ORIENT_SPECIFIC): 0 R->L, 1 L->R, 2 P->A, 3 A->P, 4 I->S, 5 S->I
    switch lower(view)
        case 'axial',    want = [0 3 4];   % across R->L, down A->P, slices counted I->S
        case 'sagittal', want = [3 5 1];   % across A->P, down S->I, slices counted L->R
        case 'coronal',  want = [0 5 3];   % across R->L, down S->I, slices counted A->P
        otherwise, error('view must be ''axial'', ''sagittal'' or ''coronal''');
    end
    perm = zeros(1,3);  flp = false(1,3);
    for q = 1:3
        a = find(floor(orientCodes(:)'/2) == floor(want(q)/2));  assert(numel(a) == 1, 'unexpected ORIENT_SPECIFIC');
        perm(q) = a;  flp(q) = orientCodes(a) ~= want(q);
    end
end
function V = toView(V, perm, flp)          % dataset array (3-D or 4-D) -> viewer layout
    V = permute(V, [perm 4]);  for q = 1:3, if flp(q), V = flip(V,q); end, end
end
function b = roiBox(roi, margin)           % [xmin xmax ymin ymax] (x = first index, y = second): ROI box over all slices + margin, whole voxels
    [bi,bj] = ind2sub([size(roi,1) size(roi,2)], find(any(roi,3)));
    b = [max(1,min(bi)-margin)-0.5  min(size(roi,1),max(bi)+margin)+0.5  max(1,min(bj)-margin)-0.5  min(size(roi,2),max(bj)+margin)+0.5];
end
function labelCells(L, nCells, box)        % cell numbers beside the strip, on the side with more room
    [si,sj] = find(L > 0);  alongI = (max(si)-min(si)) > (max(sj)-min(sj));
    for c = 1:nCells
        [ci,cj] = find(L == c);
        if alongI      % strip runs along x: numbers above or below it
            if (min(sj) - box(3)) > (box(4) - max(sj)), y = min(sj) - 1.6; else, y = max(sj) + 1.6; end
            text(mean(ci), y, num2str(c), 'Color',[0.1 0.1 0.1], 'BackgroundColor','w', 'Margin',1, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
        else           % strip runs along y: numbers left or right of it
            if (min(si) - box(1)) > (box(2) - max(si)), x = min(si) - 1.6; else, x = max(si) + 1.6; end
            text(x, mean(cj), num2str(c), 'Color',[0.1 0.1 0.1], 'BackgroundColor','w', 'Margin',1, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
        end
    end
end
function scaleBar(box, lenVox, axisDir)    % 1 mm bar in the lower left corner of the panel
    x0 = box(1) + 0.08*(box(2)-box(1));  h = 0.05*(box(4)-box(3));
    if strcmpi(axisDir,'xy'), y = box(3) + h; yt = y + 1.2*h; else, y = box(4) - h; yt = y - 1.2*h; end
    plot([x0 x0+lenVox], [y y], 'w-', 'LineWidth',4);
    text(x0 + lenVox/2, yt, '1mm', 'Color','w', 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end
