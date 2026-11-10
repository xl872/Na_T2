%% =====================================================================
%  ROI_T2long_T2vsADC.m      (2026-09-21; successor of ROI_ADChighlow_T2vsADC.m)
%  ADC vs sodium T2* (short = T2f, long = T2s) for two ROIs that are defined ON THE T2* LONG MAP and consist of
%  WHOLE SODIUM VOXELS.  Every voxel that enters the statistics therefore has ONE T2* value and
%  ONE ADC value, no sodium voxel is split between ROIs and no coverage cut of a hand-drawn proton ROI is needed.
%
%  INPUTS
%    this folder    : ROI_shrink+orig  = the DWI ROI (proton grid).  It alone defines the organoid area, for the ADC
%                                        and for the T2* maps (the same ROI masks the three maps of alignT2toDTI613New.m 6b)
%                     ADC+orig, allave+orig, ROI+orig (anatomy / zoom box of the maps)
%                     T2_to_613New_xform.mat  (alignT2toDTI613New.m step 4, including extraShift_mm)
%                     T2s_on613New+orig       (step 5a; consistency check only)
%                     ADC_R2+orig             (myADCProt.m: R^2 = gof.rsquare of the voxel-wise ADC fit, proton grid)
%    ../Na23_scan402-418 : T2allf, T2alls, T2all +orig  (T2estimate.m, whole-matrix fit, ms), allave+orig (grid info)
%                     T2allR2bi, T2allR2mono +orig   (T2estimate.m: R^2 = gof.rsquare of the voxel-wise biexponential / mono-exponential fit)
%                     ROI_T2longhigh+orig, ROI_T2longlow+orig   only for roiSource = 'file' (drawn in AFNI on T2alls)
%
%  METHOD
%    1. every DWI-ROI proton voxel is sent to its nearest sodium voxel with the stored transform (exactly how the
%       T2*_on613New maps were made).  Per sodium voxel: coverage = fraction of its ~27 proton voxels that belong to the
%       DWI ROI, ADC = mean ADC of those voxels (box average).
%       sliceOnly = true (default): ONLY the proton slice kk that is shown in the maps is used - the DWI-ROI voxels of that
%       slice, the sodium voxels seen in it, and the in-plane coverage on that slice.  The ROIs, the dots and the statistics
%       then correspond exactly to the ADC / T2* maps and the ROI map of the page (every T2* value of the plot is a block
%       that is visible in the T2* map).  sliceOnly = false = all slices (3-D), as in the versions before 2026-09-21 15:00.
%    2. eligible sodium voxels = coverage >= minCov, with a sodium fit, T2* long not on its fit bound.
%    3. ROIs from the T2* long map of the eligible voxels (roiSource = 'auto'):
%          splitMode 'median'  : T2* long above / below the median (all eligible voxels are used)
%                    'tertile' : upper / lower third (middle third unused)
%                    [pLow pHigh]: <= pLow-th / >= pHigh-th percentile
%       or centre + ring (roiSource = 'centre'): T2* long HIGH = the central highest T2* long voxel of the sodium slice that
%       is shown in the maps (eligible voxel with the highest T2* long within centreRad of the DWI-ROI centroid, or centreVox),
%       T2* long LOW = its 8 in-plane neighbours (the ring around it); ring voxels only need a DWI-ROI coverage of
%       ringMinCov (rim voxels are partly outside the DWI ROI), so that the whole ring is used;
%       or hand-drawn sodium-grid ROIs (roiSource = 'file'); their voxels must be eligible.
%       The ADC is NOT used to choose the voxels.  So the T2* long difference between the ROIs is by construction
%       (no test is reported for it); the ADC, T2* short and T2* mono differences are the results.
%    4. per ROI and quantity: mean, SD, SE = SD/sqrt(n); two-sample t-test high vs low.
%       T2* (short, long, mono): one value per SODIUM voxel, n = sodium voxels.
%       ADC with adcUnit = 'proton' (default): one value per ADC-grid (proton) voxel of the DWI ROI inside the ROI's sodium
%       voxels, n = proton voxels - so a ROI of one sodium voxel still has ~27 ADC values, a mean and an error bar.  The plot
%       shows every ADC-grid voxel as a transparent dot (x = its own ADC, y = T2* of the sodium voxel it lies in) and the
%       mean +/- SE on top.  Neighbouring ADC voxels are NOT independent (smoothed DWI, same tissue), so the per-voxel ADC
%       p value is too small; 'ADCbox' (one box-averaged ADC per sodium voxel, as before) is kept in the tables for reference.
%       A ROI of ONE voxel has no SD / SE (no error bar); the pooled-variance t-test is then the single-case test of
%       Crawford & Howell: t = (x - mean_ring) / (SD_ring * sqrt(1 + 1/n)), df = n - 1.
%       NOTE: extraShift_mm of the alignment was chosen so that the central T2* voxel lies on the high-ADC area, so the
%       ADC of the centre voxel is not independent of that choice (the ADC p value is optimistic).
%       The T2* short values drop the voxels whose T2f sits on its own fit bound (dropOwnBound), as before.
%    5. fit quality (the "fitting R2" rows of the supplementary tables): the R^2 (gof.rsquare) of the voxel-wise ADC fit
%       (ADC_R2, per ADC voxel = the same voxels as the ADC row) and of the voxel-wise biexponential and mono-exponential
%       T2* fits (T2*R2bi / T2*R2mono, per sodium voxel; ALL voxels of the ROI, no bound filter) are summarised like the
%       other quantities (mean, SD, SE, n) and written to the stats CSV as ADC_R2, T2bi_R2 and T2mono_R2.
%
%  OUTPUTS
%    this folder    : ROI_T2long_T2vsADC.jpg / .eps        T2* short vs ADC, T2* long vs ADC (mean +/- SE per ROI)
%                     ROI_T2long_map_protongrid.jpg        the two ROIs as solid fills on the proton anatomy (slice kk)
%                     ROI_T2long_T2vsADC_stats.csv / _ttest.csv,  ROI_T2long_sodiumvoxels.csv (all eligible voxels)
%                        (the stats CSV also holds the fit-quality rows ADC_R2, T2bi_R2, T2mono_R2; the voxel CSV the R^2 per voxel)
%                     ROI_T2longhigh+orig, ROI_T2longlow+orig    the ROIs on the proton grid (DWI-ROI voxels of those sodium voxels)
%    ../Na23_scan402-418 : ROIlabel_T2long_on402grid+orig       1 = T2* long high, 2 = low, 3 = eligible but unused
%
%  Run after alignT2toDTI613New.m (the script changes to data/DWI_scan613 itself), cell by cell (Ctrl+Enter) or as a whole.
%  =====================================================================
cd(fullfile(fileparts(mfilename('fullpath')),'..','data','DWI_scan613'));   % release layout: run in the dataset folder
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','helper'));   % release layout: readT2Bounds.m
clc; clear; close all;

%% ---------------- 0. settings ----------------------------------------------------
here = pwd;  T2D = fullfile(fileparts(pwd),'Na23_scan402-418');
orgName   = 'Fig. 5 #1, scan 613';   orgShort = 'Fig. 5 #1 (613)';   % console / figure title
xformFile = 'T2_to_613New_xform.mat';   chkFile = 'T2s_on613New';   gridTag = '402grid';
dwiROI    = 'ROI_shrink';   % the DWI ROI in this folder (proton grid); the same one alignT2toDTI613New.m uses for its maps
T2names   = {'T2allf','T2alls','T2all'};   whichFit = 'last';   % fast, slow, mono maps in T2D + their fit in T2estimate.m
R2names   = {'T2allR2bi','T2allR2mono'};   % R^2 maps of those fits in T2D (T2estimate.m, gof.rsquare): biexponential, mono-exponential
adcR2name = 'ADC_R2';       % R^2 map of the voxel-wise ADC fit in this folder (myADCProt.m, gof.rsquare)
roiSource = 'centre';       % 'centre' = central highest T2* long voxel vs the ring of its in-plane neighbours (see METHOD 3);
                            % 'auto' = split the eligible sodium voxels by T2* long (splitMode);  'file' = roiFilesT
splitMode = 'median';       % roiSource 'auto': 'median' | 'tertile' | [pLow pHigh] (percentiles of T2* long over the eligible voxels)
centreSlices = [];          % roiSource 'centre': sodium slices (k, 1-based) to use; [] = the sodium slice shown in proton slice kk
centreVox = [];             %   centre voxel per slice as rows [i j k] (1-based sodium indices); [] = automatic (see centreRad)
ringMinCov = 0.25;          %   ring voxels: DWI-ROI coverage needed (lower than minCov: the corner voxels of the ring reach the rim,
                            %   e.g. the bottom-left one of scan 613 is covered 0.33); only their DWI-ROI part enters the ADC
centreRad = 1.0;            %   automatic centre = highest T2* long among the eligible voxels within this in-plane distance
                            %   (sodium voxels) of the DWI-ROI centroid
roiFilesT = {'ROI_T2longhigh','ROI_T2longlow'};   % sodium-grid ROIs in T2D (read for roiSource = 'file'); also the names
                                                  % of the proton-grid ROI datasets written into this folder
roiLabels = {'T2* long high','T2* long low'};
minCov    = 0.5;            % sodium voxel belongs to the organoid when the DWI ROI covers at least this fraction of it
dropOwnBound = true;        % T2* short: drop sodium voxels whose T2f sits on its fit bound (T2* long bound voxels are not eligible)
T2fBounds = [];  T2sBounds = [];   % [] = read from T2estimate_scan402_418.m with readT2Bounds.m; or e.g. [0.1 2] / [2 50]
boundTol  = 1e-6;           % "on a bound" = within this distance (ms) of a bound
vartype   = 'equal';        % t-test: 'equal' = ordinary two-sample (Student) t-test (default); 'unequal' = Welch
adcUnit   = 'proton';       % 'proton' = ADC statistics and dots per ADC-grid (proton) voxel; 'sodium' = per sodium voxel (box average)
showVoxels = true;          % draw the individual voxels (adcUnit) as transparent dots behind the mean +/- SE
dotSize = 30;  dotAlpha = 0.30;     % size / opacity of those dots
errBar    = 'SE';           % error bars of the mean: 'SE' (= SD/sqrt(n), as on the earlier pages) or 'SD'
kk        = 14;             % coronal proton slice of the ROI map (the slice of alignT2toDTI613New.m 6b)
sliceOnly = true;           % true = analyse ONLY proton slice kk (ROIs, dots and statistics match the maps shown); false = all slices
mapMargin = 3;              % ROI map zoom: organoid ROI box + this many proton voxels, as in alignT2toDTI613New.m 6b; [] = full FOV
mapAlpha  = 1;              % opacity of the ROI colour fill (1 = solid)
col = [0.165 0.471 0.839;   % T2* long high : blue
       0.922 0.408 0.204];  % T2* long low  : orange
ink = [0 0 0];              % text / axis colour (black, as in Final Figures.pptx)
fs  = 20;  axLW = 2;        % plot font size (pt, Arial bold) and axis line width
legLabels = {'ROI1','ROI2'};% plot legend (blue / orange dot) for T2* long high / low
showTitles = false;         % true: panel titles with the t-test p values (they were cropped away on the slide)
adcScale = 100;             % ADC axis in 10^-4 mm^2/s (= um^2/s / 100)
yLimShort = [0 1.2];  yTickShort = 0:0.4:1.2;   % T2* short axis; [] = automatic
yLimLong  = [3 10];   yTickLong  = 4:2:10;      % T2* long axis;  [] = automatic
nTickMax = 4;               % at most this many y ticks on the automatic axes

%% ---------------- 1. load --------------------------------------------------------
[adc,InfoD] = BrikLoad(fullfile(here,'ADC+orig.BRIK'));  adc = double(adc);  ND = size(adc);
anatD = double(BrikLoad(fullfile(here,'allave+orig.BRIK')));
roiBox = BrikLoad(fullfile(here,'ROI+orig.BRIK')) > 0;                   % zoom box of the maps only
roiD  = BrikLoad(fullfile(here,[dwiROI '+orig.BRIK'])) > 0;              % the DWI ROI
T2f  = double(BrikLoad(fullfile(T2D,[T2names{1} '+orig.BRIK'])));
T2s  = double(BrikLoad(fullfile(T2D,[T2names{2} '+orig.BRIK'])));
T2m  = double(BrikLoad(fullfile(T2D,[T2names{3} '+orig.BRIK'])));
R2bi = double(BrikLoad(fullfile(T2D,[R2names{1} '+orig.BRIK'])));    % R^2 of the biexponential fit per sodium voxel (T2estimate.m)
R2mo = double(BrikLoad(fullfile(T2D,[R2names{2} '+orig.BRIK'])));    % R^2 of the mono-exponential fit
adcR2 = double(BrikLoad(fullfile(here,[adcR2name '+orig.BRIK'])));   % R^2 of the ADC fit per proton voxel (myADCProt.m)
[~,InfoT] = BrikLoad(fullfile(T2D,'allave+orig.BRIK'));  NT = size(T2f);
X = load(fullfile(here,xformFile),'M_Dijk_to_T2ijk','shiftMode','extraShift_mm');
M = X.M_Dijk_to_T2ijk;                                  % proton ijk (0-based) -> sodium ijk (0-based)
resD = abs(InfoD.DELTA(1));  resT = abs(InfoT.DELTA(1));
fprintf('%s\nproton grid %dx%dx%d @ %.2f mm, sodium grid %dx%dx%d @ %.2f mm; transform: shiftMode=%s extraShift=[%g %g %g] mm\n', ...
        orgName, ND, resD, NT, resT, X.shiftMode, X.extraShift_mm);
fprintf('DWI ROI %s: %d proton voxels (%.2f mm^3 = %.1f sodium-voxel volumes), %d without an ADC fit\n', ...
        dwiROI, nnz(roiD), nnz(roiD)*resD^3, nnz(roiD)*resD^3/resT^3, nnz(roiD & adc<=0));
if isempty(T2fBounds) || isempty(T2sBounds)
    [T2fBounds, T2sBounds, bnd] = readT2Bounds(fullfile('..','..','scripts','T2estimate_scan402_418.m'), whichFit);
    fprintf('fit bounds read from T2estimate.m line %d (%s): T2f [%g %g], T2s [%g %g] ms\n', bnd.line, bnd.model, T2fBounds, T2sBounds);
end
% the R^2 maps must come from the same fits as the parameter maps: same voxels, and not older than them
if any(T2s(:) > 0 & R2bi(:) == 0) || any(adc(:) > 0 & adcR2(:) == 0)
    warning('R^2 maps and parameter maps do not cover the same voxels - rerun T2estimate.m / myADCProt.m so that both come from one fit');
end
dR = dir(fullfile(T2D,[R2names{1} '+orig.HEAD']));  dT = dir(fullfile(T2D,[T2names{2} '+orig.HEAD']));
if dR.datenum < dT.datenum, warning('%s is older than %s: the R^2 map may not belong to the current T2* maps', R2names{1}, T2names{2}); end
dR = dir(fullfile(here,[adcR2name '+orig.HEAD']));  dA = dir(fullfile(here,'ADC+orig.HEAD'));
if dR.datenum < dA.datenum, warning('%s is older than ADC: the R^2 map may not belong to the current ADC map', adcR2name); end
onBf = abs(T2f-T2fBounds(1))<boundTol | abs(T2f-T2fBounds(2))<boundTol;    % T2f on its own bound
onBs = abs(T2s-T2sBounds(1))<boundTol | abs(T2s-T2sBounds(2))<boundTol;    % T2s on its own bound
hasFit = T2s(:) > 0 & T2m(:) > 0;                                           % sodium voxels with a fit

%% ---------------- 2. DWI ROI and ADC onto the sodium grid -------------------------------
%  every proton voxel -> nearest sodium voxel (same transform that made T2*_on613New)
[I,J,K] = ndgrid(0:ND(1)-1, 0:ND(2)-1, 0:ND(3)-1);
tq = M*[I(:)'; J(:)'; K(:)'; ones(1,numel(I))];  tr = round(tq(1:3,:));  clear I J K tq
in = (all(tr>=0,1) & tr(1,:)<=NT(1)-1 & tr(2,:)<=NT(2)-1 & tr(3,:)<=NT(3)-1)';
lin = zeros(numel(in),1);  lin(in) = sub2ind(NT, tr(1,in)+1, tr(2,in)+1, tr(3,in)+1);  clear tr
inS = in;                                                            % proton voxels that take part
if sliceOnly, [~,~,kpAll] = ind2sub(ND,(1:prod(ND))');  inS = in & kpAll == kk;  clear kpAll; end      % ... only those of the shown slice
cntAll = accumarray(lin(inS), 1, [prod(NT) 1]);                      % proton voxels per sodium voxel (~27)
chk = double(BrikLoad(fullfile(here,[chkFile '+orig.BRIK'])));      % consistency with the maps of the alignment script
d = max(abs(chk(in) - T2s(lin(in))));
if d > 1e-3, warning('%s differs from the stored transform / %s (max %.3g ms) - rerun the alignment script, steps 4-5', chkFile, T2names{2}, d); end

m   = roiD(:) & inS & adc(:) > 0;                                     % DWI-ROI proton voxels with an ADC fit
cnt = accumarray(lin(m), 1,      [prod(NT) 1]);
cov = cnt ./ max(cntAll,1);                                           % coverage of each sodium voxel by the DWI ROI
adcBox = accumarray(lin(m), adc(m), [prod(NT) 1]) ./ max(cnt,1);      % box-averaged ADC of the DWI-ROI voxels in it
elig = cov >= minCov & hasFit & ~onBs(:);
fprintf('%s: sodium voxels touched by the DWI ROI: %d; coverage >= %.2f: %d; without a fit %d, T2* long on its bound %d -> %d eligible\n', ...
        ternary(sliceOnly,sprintf('proton slice k=%d only',kk),'all slices'), nnz(cov>0), minCov, nnz(cov>=minCov), nnz(cov>=minCov & ~hasFit), nnz(cov>=minCov & hasFit & onBs(:)), nnz(elig));
for c = [0.35 0.5 0.65 0.8 0.95], fprintf('   coverage >= %.2f : %2d voxels\n', c, nnz(cov>=c & hasFit & ~onBs(:))); end
E = find(elig);
fprintf('Spearman over the %d eligible sodium voxels, ADC vs:', numel(E));
for q = {T2s,'T2* long'; T2f,'T2* short'; T2m,'T2* mono'}'
    [rs,ps] = corr(adcBox(E), q{1}(E), 'type','Spearman');  fprintf('   %s rho = %+.2f (p = %.3f)', q{2}, rs, ps);
end
fprintf('\n');

%% ---------------- 3. the two ROIs from the T2* long map (whole sodium voxels) --------------
assign = zeros(prod(NT),1);                       % 1 = T2* long high, 2 = low, 3 = eligible, unused
switch roiSource
    case 'auto'
        v = T2s(E);
        if ischar(splitMode) || isstring(splitMode)
            switch lower(char(splitMode))
                case 'median',  cutLo = median(v);  cutHi = cutLo;  hi = v > cutHi;  lo = v < cutLo;
                case 'tertile', cu = quantile(v,[1/3 2/3]);  cutLo = cu(1); cutHi = cu(2);  hi = v > cutHi;  lo = v <= cutLo;
                otherwise, error('splitMode must be median, tertile or [pLow pHigh]');
            end
            ruleTxt = char(splitMode);
        else
            cutLo = prctile(v,splitMode(1));  cutHi = prctile(v,splitMode(2));  hi = v >= cutHi;  lo = v <= cutLo;
            ruleTxt = sprintf('P%g / P%g', splitMode);
        end
        assign(E) = 3;  assign(E(hi)) = 1;  assign(E(lo)) = 2;
        fprintf('ROIs from the T2* long map (%s split of the eligible voxels): high > %.2f ms, low <%s %.2f ms\n', ruleTxt, cutHi, ternary(cutLo<cutHi,'=',''), cutLo);
    case 'centre'
        assign(E) = 3;
        [gi,gj,gk] = ndgrid(1:NT(1),1:NT(2),1:NT(3));
        cen = [sum(cov.*gi(:)) sum(cov.*gj(:)) sum(cov.*gk(:))] / sum(cov);        % DWI-ROI centroid on the sodium grid
        if isempty(centreSlices)                                                   % sodium slice seen in proton slice kk
            [~,~,kp] = ind2sub(ND, find(m));  [~,~,ks] = ind2sub(NT, lin(m));  centreSlices = mode(ks(kp == kk));
        end
        fprintf('centre + ring ROIs: DWI-ROI centroid on the sodium grid [%.2f %.2f %.2f], sodium slice(s) %s\n', cen, mat2str(centreSlices));
        eligT = reshape(elig,NT);  cvox = zeros(0,3);
        ringOK = reshape(cov >= ringMinCov & hasFit & ~onBs(:), NT);               % looser coverage cut for the ring
        for sl = centreSlices(:)'
            if ~isempty(centreVox), c = centreVox(centreVox(:,3) == sl,:);
            else
                cand = find(eligT & gk == sl & hypot(gi-cen(1), gj-cen(2)) <= centreRad);
                [~,o] = max(T2s(cand));  [c1,c2,c3] = ind2sub(NT, cand(o));  c = [c1 c2 c3];
            end
            if isempty(c), warning('no centre voxel in sodium slice %d', sl); else, cvox = [cvox; c]; end %#ok<AGROW>
        end
        cl = sub2ind(NT, cvox(:,1), cvox(:,2), cvox(:,3));
        for c = 1:size(cvox,1)                                                     % ring = the 8 in-plane neighbours
            ring = ringOK & gk == cvox(c,3) & max(abs(gi-cvox(c,1)), abs(gj-cvox(c,2))) == 1;
            ring(cl) = false;  assign(ring(:) & assign ~= 1) = 2;  assign(cl(c)) = 1;
            nb = gk == cvox(c,3) & max(abs(gi-cvox(c,1)), abs(gj-cvox(c,2))) == 1;
            fprintf('   centre [%d %d %d]: T2* long %.1f ms, coverage %.2f; ring: %d of 8 neighbours used (coverage >= %.2f: %.2f - %.2f; T2* long %.1f - %.1f ms)\n', ...
                    cvox(c,:), T2s(cl(c)), cov(cl(c)), nnz(ring), ringMinCov, min(cov(ring)), max(cov(ring)), min(T2s(ring)), max(T2s(ring)));
            for q = find(nb & ~ringOK)', [q1,q2,q3] = ind2sub(NT,q); fprintf('      neighbour [%d %d %d] not used (coverage %.2f)\n', q1, q2, q3, cov(q)); end
        end
        ruleTxt = 'centre + ring';  clear gi gj gk
    case 'file'
        assign(E) = 3;
        for r = 1:2
            rT = BrikLoad(fullfile(T2D,[roiFilesT{r} '+orig.BRIK'])) > 0;
            fprintf('%s: %d sodium voxels drawn, %d of them eligible (the others are dropped)\n', roiFilesT{r}, nnz(rT), nnz(rT(:) & elig));
            if any(assign(rT(:) & elig) < 3), warning('%s overlaps %s; shared voxels stay with %s', roiFilesT{r}, roiFilesT{1}, roiFilesT{1}); end
            assign(rT(:) & elig & assign == 3) = r;
        end
        ruleTxt = 'hand-drawn';
    otherwise, error('roiSource must be centre, auto or file');
end
nR = 2;  S = struct('label',{},'lin',{});
for r = 1:nR
    S(r).label = roiLabels{r};  S(r).lin = find(assign == r);
    fprintf('   %-14s %2d sodium voxels  (%s; T2* long %.1f - %.1f ms, mean coverage %.2f)\n', roiLabels{r}, numel(S(r).lin), ruleTxt, min(T2s(S(r).lin)), max(T2s(S(r).lin)), mean(cov(S(r).lin)));
end
writeOut(reshape(assign,NT), InfoT, fullfile(T2D,['ROIlabel_T2long_on' gridTag]));
roiP = cell(1,nR);                                  % the ROIs on the proton grid = DWI-ROI voxels of those sodium voxels
for r = 1:nR
    a = false(ND);  a(m) = assign(lin(m)) == r;  roiP{r} = a;
    writeOut(double(a), InfoD, fullfile(here,roiFilesT{r}));
end

%% ---------------- 4. statistics ----------------------------------------------------------
%  T2* : one value per sodium voxel.   ADC : per ADC-grid (proton) voxel (adcUnit = 'proton') or per sodium voxel (box average).
%  ADCbox = the box-averaged ADC per sodium voxel, always listed for reference.
%  Fit quality (R2 rows of the tables): ADC_R2 per ADC voxel (same voxels as ADC), T2bi_R2 / T2mono_R2 per sodium voxel (all ROI voxels).
qty  = {'ADC','T2short','T2long','T2mono','ADCbox','ADC_R2','T2bi_R2','T2mono_R2'};
qlab = {'ADC (\mum^2/s)','T2* short component (ms)','T2* long component (ms)','T2* mono (ms)','ADC box average (\mum^2/s)', ...
        'ADC fit R^2','T2* biexponential fit R^2','T2* mono-exponential fit R^2'};
adcR2Box = accumarray(lin(m), adcR2(m), [prod(NT) 1]) ./ max(cnt,1);          % box-averaged ADC-fit R^2 per sodium voxel (adcUnit = 'sodium')
maps = {reshape(adcBox,NT),T2f,T2s,T2m,reshape(adcBox,NT),reshape(adcR2Box,NT),R2bi,R2mo};
onB  = repmat({false(NT)},1,numel(qty));  onB{2} = onBf;                       % only T2short drops the voxels on the T2f bound
iADC_R2 = 6;
perProton = strcmp(adcUnit,'proton');  unitQ = repmat({'sodium voxel'},1,numel(qty));  if perProton, unitQ([1 iADC_R2]) = {'ADC voxel'}; end
nQ = numel(qty);  mn = nan(nR,nQ); sd = mn; se = mn; nUse = zeros(nR,nQ); nDrop = nUse;  vals = cell(nR,nQ);
Xv = cell(nR,nQ);  Yv = cell(nR,nQ);  mnX = nan(nR,nQ);  seX = mnX;  sdX = mnX;      % per panel: dots (x = ADC, y = T2*) and the ADC mean / SE of the SAME voxels
for r = 1:nR
    for q = 1:nQ
        keep = true(size(S(r).lin));
        if dropOwnBound, keep = ~onB{q}(S(r).lin); end
        L = S(r).lin(keep);  v = maps{q}(L);
        if perProton, pm = m & ismember(lin, L);  x = adc(pm);  y = maps{q}(lin(pm));      % the DWI-ROI proton voxels inside these sodium voxels
        else,         x = adcBox(L);  y = maps{q}(L);  end
        if q == 1, v = x; end                                                              % ADC in its own unit
        if q == iADC_R2 && perProton, v = adcR2(pm); end                                   % ADC-fit R^2 per ADC voxel (the voxels of the ADC row)
        Xv{r,q} = x(:);  Yv{r,q} = y(:);  mnX(r,q) = mean(x);  if numel(x) >= 2, sdX(r,q) = std(x);  seX(r,q) = sdX(r,q)/sqrt(numel(x)); end
        vals{r,q} = v(:);  nUse(r,q) = numel(v);  nDrop(r,q) = nnz(~keep);
        if ~isempty(v), mn(r,q) = mean(v); end
        if numel(v) >= 2, sd(r,q) = std(v);  se(r,q) = sd(r,q)/sqrt(numel(v)); end      % one voxel: no SD / SE
    end
end
fprintf('\nmean +/- SE (SD, n); ADC per %s, T2* per sodium voxel%s\n', unitQ{1}, ternary(dropOwnBound,'; T2short without the voxels on the T2f fit bound',''));
fprintf('%-16s', 'ROI'); fprintf('%33s', qty{:}); fprintf('\n');
for r = 1:nR
    fprintf('%-16s', roiLabels{r});
    for q = 1:nQ, fprintf('   %8.2f +/- %6.2f (%6.2f, n=%3d)', mn(r,q), se(r,q), sd(r,q), nUse(r,q)); end
    fprintf('\n');
end
fprintf('  T2short: sodium voxels dropped for sitting on the T2f fit bound  %s %d, %s %d\n', roiLabels{1}, nDrop(1,2), roiLabels{2}, nDrop(2,2));
fprintf('  fit quality for the R2 rows of the table (mean +/- SD, n): ADC fit R^2 per %s, T2* fit R^2 per sodium voxel (all ROI voxels)\n', unitQ{iADC_R2});
for r = 1:nR
    fprintf('    %-14s ADC R2 %.3f +/- %.3f (n=%d)   T2 biexp R2 %.3f +/- %.3f (n=%d)   T2 mono R2 %.3f +/- %.3f (n=%d)\n', roiLabels{r}, ...
            mn(r,6), sd(r,6), nUse(r,6), mn(r,7), sd(r,7), nUse(r,7), mn(r,8), sd(r,8), nUse(r,8));
end

[rr,qq] = ndgrid(1:nR,1:nQ);
Tstat = table(roiLabels(rr(:))', qty(qq(:))', unitQ(qq(:))', mn(:), sd(:), se(:), nUse(:), nDrop(:), ...
    'VariableNames',{'ROI','quantity','unit','mean','SD','SE','n','n_sodiumVoxelsDroppedOnBound'});
writetable(sortrows(Tstat,{'ROI','quantity'}), fullfile(here,'ROI_T2long_T2vsADC_stats.csv'));
Ea = find(assign > 0);  [si,sj,sk] = ind2sub(NT, Ea);  lab = {'high','low','unused'};       % eligible voxels + the ring voxels
Tv = table(lab(assign(Ea))', si, sj, sk, cov(Ea), cnt(Ea), adcBox(Ea), T2f(Ea), T2s(Ea), T2m(Ea), onBf(Ea), R2bi(Ea), R2mo(Ea), adcR2Box(Ea), ...
     'VariableNames',{'ROI_T2long','i','j','k','DWIROIcoverage','n_ADCvoxels','ADCbox_1e-6mm2s','T2fast_ms','T2slow_ms','T2mono_ms','T2f_onBound', ...
                      'T2bi_R2','T2mono_R2','ADC_R2box'});
writetable(sortrows(Tv,'T2slow_ms','descend'), fullfile(here,'ROI_T2long_sodiumvoxels.csv'));
fprintf('  wrote ROI_T2long_T2vsADC_stats.csv and ROI_T2long_sodiumvoxels.csv\n');

%% ---------------- 4b. two-sample t-test: T2* long high vs low -----------------------------
tname = ternary(strcmp(vartype,'unequal'),'Welch t-test','two-sample t-test');
pT = nan(1,nQ); single = false; Tt = table('Size',[nQ 10], 'VariableTypes',[{'string'} repmat({'double'},1,9)], ...
     'VariableNames',{'quantity','n_high','n_low','diff_high_minus_low','t','df','p_used','p_welch','p_student','p_ranksum'});
fprintf('\n%s %s vs %s; ADC per %s, the others per sodium voxel (T2long: the ROIs were chosen on this map, its p value is by construction):\n', tname, roiLabels{1}, roiLabels{2}, unitQ{1});
fprintf('  %-8s %5s %5s %10s %8s %7s %10s %10s %12s\n', 'quantity','n_hi','n_lo','diff','t','df','p_Student','p_Welch','p_MannWhitney');
for q = 1:nQ
    a = vals{1,q}; b = vals{2,q}; t = NaN; df = NaN; pW = NaN; pS = NaN; pR = NaN;
    if numel(a) >= 2 && numel(b) >= 2
        [~,pW,~,sW] = ttest2(a,b,'Vartype','unequal');  [~,pS,~,sS] = ttest2(a,b,'Vartype','equal');  pR = ranksum(a,b);
        st = ternary(strcmp(vartype,'unequal'),sW,sS);  t = st.tstat; df = st.df;
        pT(q) = ternary(strcmp(vartype,'unequal'),pW,pS);
    elseif min(numel(a),numel(b)) == 1 && max(numel(a),numel(b)) >= 3
        [~,pS,~,sS] = ttest2(a,b,'Vartype','equal');  t = sS.tstat; df = sS.df;  pR = ranksum(a,b);   % single voxel vs a group:
        pT(q) = pS;  single = true;                        % pooled-variance t-test = Crawford-Howell single-case test (no Welch)
    end
    Tt(q,:) = {qty{q}, numel(a), numel(b), mean(a)-mean(b), t, df, pT(q), pW, pS, pR};
    fprintf('  %-8s %5d %5d %10.2f %8.2f %7.1f %10.4f %10.4f %12.4f%s\n', qty{q}, numel(a), numel(b), mean(a)-mean(b), t, df, pS, pW, pR, ...
            ternary(q==3,'   (by construction)', ternary((q==1 || q==iADC_R2) && perProton,'   (ADC voxels are not independent: p too small)','')));
end
if single, tname = sprintf('t-test (%d vs %d voxels)', numel(S(1).lin), numel(S(2).lin)); fprintf('  one ROI is a single voxel: Student t = single-case test (Crawford & Howell), SD of the other ROI, df = n - 1\n'); end
writetable(Tt, fullfile(here,'ROI_T2long_T2vsADC_ttest.csv'));
fprintf('  wrote ROI_T2long_T2vsADC_ttest.csv\n');

%% ---------------- 5. plot: T2* short vs ADC, T2* long vs ADC ------------------------
figure('Name','ROIs from the T2* long map: T2* vs ADC (sodium grid)','Position',[100 100 924 400],'Color','w');   % aspect of the slide frame (5.79 x 2.51 in)
panels = [2 3];                                             % T2short, T2long
if strcmpi(errBar,'SD'), eY = sd; eX = sdX; else, eY = se; eX = seX; end      % error bars: SD or SE
xv = [mnX(:,panels)-eX(:,panels); mnX(:,panels)+eX(:,panels); mnX(:,panels)];  xv = xv(:);
if showVoxels, xv = [xv; vertcat(Xv{:,panels})]; end
xv = xv(~isnan(xv));  xl = [min(xv) max(xv)];  xl = xl + [-1 1]*max(0.08*diff(xl), 5);   % common ADC axis
pLab = {'','T_2*_{short} (ms)','T_2*_{long} (ms)'};     % y labels written as in Final Figures.pptx
for p = 1:2
    q = panels(p); ax = subplot(1,2,p); hold(ax,'on');
    h = gobjects(1,nR);
    if showVoxels                                           % every voxel (adcUnit): its own ADC vs the T2* of its sodium voxel
        for r = 1:nR, scatter(Xv{r,q}/adcScale, Yv{r,q}, dotSize, col(r,:), 'filled', 'MarkerFaceAlpha',dotAlpha, 'MarkerEdgeColor','none'); end
    end
    for r = 1:nR                                            % mean +/- SE on top: x over the ADC values, y over the sodium voxels
        e = [eY(r,q) eX(r,q)/adcScale];  e(isnan(e)) = 0;  % single-voxel ROI: no T2* error bar
        errorbar(mnX(r,q)/adcScale, mn(r,q), e(1), e(1), e(2), e(2), 'o', ...
                        'Color',col(r,:), 'MarkerFaceColor',col(r,:), 'MarkerEdgeColor','w', ...
                        'MarkerSize',9, 'LineWidth',1.5, 'CapSize',ternary(any(e),8,0));
    end
    yv = [mn(:,q); mn(:,q)-eY(:,q); mn(:,q)+eY(:,q)];  if showVoxels, yv = [yv; vertcat(Yv{:,q})]; end
    yv = yv(~isnan(yv));
    if ~isempty(yv), yl = [min(yv) max(yv)]; yl = yl + [-1 1]*max(0.15*diff(yl), 0.05*abs(mean(yl))); ylim(yl); end
    if q == 2 && ~isempty(yLimShort), ylim(yLimShort); yticks(yTickShort);
    elseif q == 3 && ~isempty(yLimLong), ylim(yLimLong); yticks(yTickLong);
    else, yticks(fewTicks(ylim, nTickMax)); end
    xlim(xl/adcScale);  xticks(fewTicks(xl/adcScale, nTickMax));
    xlabel('ADC (10^{-4} mm^2/s)'); ylabel(pLab{q});
    if showTitles
        if q == 3, l2 = 'defines the ROIs (no test)'; else, l2 = sprintf('%s: %s', tname, pstr(pT(q))); end
        title({sprintf('%s vs ADC', qlab{q}), l2}, 'FontSize',10, 'FontWeight','normal');
    end
    for r = 1:nR, h(r) = plot(NaN, NaN, 'o', 'MarkerFaceColor',col(r,:), 'MarkerEdgeColor','none', 'MarkerSize',10); end   % legend: plain dots
    legend(h, legLabels, 'Location','northoutside', 'Orientation','horizontal', 'Box','off', 'TextColor',ink);
    styleAx(ax, ink, axLW);
end
if showTitles, sgtitle(sprintf('%s: mean \\pm %s; ADC n = %d / %d, %s', orgShort, upper(errBar), nUse(1,1), nUse(2,1), pstr(pT(1))), 'FontSize',11, 'Color',ink); end
fontsize(gcf,fs,"points")
set(findall(gcf,'-property','FontName'),'FontName','Arial');  set(findall(gcf,'-property','FontWeight'),'FontWeight','bold');
set(findall(gcf,'Type','axes'),'Toolbar',[]);              % no axes-toolbar icons in a batch-mode export
print(gcf, fullfile(here,'ROI_T2long_T2vsADC.jpg'), '-djpeg', '-r300');
print(gcf, fullfile(here,'ROI_T2long_T2vsADC.eps'), '-depsc', '-r300');
fprintf('  wrote ROI_T2long_T2vsADC.jpg / .eps\n');

%% ---------------- 6. map: where the two ROIs are (proton grid, slice kk) ----------------------
%  grey proton anatomy with each ROI as a solid colour fill (blue = T2* long high, orange = low), no colour bar.
%  The fills are whole sodium voxels (3x3 proton voxels in plane) clipped by the DWI ROI.
[bi,bj] = ind2sub(ND(1:2), find(any(roiBox,3)));
figure('Name','ROI maps (proton grid)','Position',[50 50 600 300],'Color','w');
for r = 1:nR
    subplot(1,2,r); showAnat(anatD(:,:,kk)); overlayMask(roiP{r}(:,:,kk), col(r,:), mapAlpha);
    if ~isempty(mapMargin), xlim([max(1,min(bi)-mapMargin)-0.5 min(ND(1),max(bi)+mapMargin)+0.5]); ylim([max(1,min(bj)-mapMargin)-0.5 min(ND(2),max(bj)+mapMargin)+0.5]); end   % whole voxels
    title(sprintf('%s, coronal k=%d', roiLabels{r}, kk), 'FontSize',10, 'FontWeight','normal');
end
sgtitle(sprintf('proton grid (%.2g mm): blue = %s, orange = %s', resD, roiLabels{1}, roiLabels{2}), 'FontSize',10, 'Color',ink);
set(findall(gcf,'Type','axes'),'Toolbar',[]);
print(gcf, fullfile(here,'ROI_T2long_map_protongrid.jpg'), '-djpeg', '-r300');
fprintf('  wrote ROI_T2long_map_protongrid.jpg (slice k=%d: %d / %d proton voxels)\n', kk, nnz(roiP{1}(:,:,kk)), nnz(roiP{2}(:,:,kk)));


%% ================= local functions ===================================================
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
function styleAx(ax, ink, lw)              % plot style of Final Figures.pptx: no grid, left + bottom axes only, thick black axes, black bold Arial
    grid(ax,'off');  box(ax,'off');
    set(ax, 'XColor',ink, 'YColor',ink, 'LineWidth',lw, 'FontName','Arial', 'FontWeight','bold');
    set([ax.XLabel ax.YLabel ax.Title], 'Color',ink, 'FontName','Arial', 'FontWeight','bold');
end
function t = fewTicks(lim, nMax)           % round tick values inside lim, at most nMax of them
    st = sort(reshape([1 2 2.5 5]'*10.^(-3:4), 1, []));
    for s = st, t = ceil(lim(1)/s-1e-9)*s : s : floor(lim(2)/s+1e-9)*s;  if numel(t) <= nMax, return; end, end
end
function s = pstr(p), if isnan(p), s = 'p = n/a'; elseif p < 0.001, s = 'p < 0.001'; else, s = sprintf('p = %.3f',p); end, end
function writeOut(Mp, Info, prefixPath)   % write a float dataset on the grid of Info (as in the alignment script)
    [p,f] = fileparts(prefixPath); old = dir(fullfile(p,[f '+orig.*'])); for q = 1:numel(old), delete(fullfile(p,old(q).name)); end
    Info.BRICK_TYPES = 3; Info.DATASET_RANK(2) = 1; Info.BRICK_FLOAT_FACS = 0;
    Info.BRICK_STATS = [min(Mp(:)) max(Mp(:))];
    for fld = {'TAXIS_NUMS','TAXIS_FLOATS','TAXIS_OFFSETS','BRICK_STATAUX','BRICK_KEYWORDS'}
        if isfield(Info,fld{1}), Info = rmfield(Info,fld{1}); end
    end
    if isfield(Info,'BRICK_LABS'), Info.BRICK_LABS = f; end
    Opt.Prefix = fullfile(p,f); Opt.OverWrite = 'y'; Opt.verbose = 0; Opt.Scale = 0; Opt.View = '+orig';
    [err,ErrMessage] = WriteBrik(double(Mp),Info,Opt); if err, error(ErrMessage); end
    fprintf('  wrote %s+orig (in %s)\n', f, p);
end
function showAnat(anat)                   % grey anatomy exactly as showMap of the alignment script draws it
    B = repmat(anat'/prctile(anat(:),99.5),1,1,3); B(B>1) = 1; image(B); axis image xy off; hold on
end
function overlayMask(m, rgb, alpha)       % filled overlay of a 2-D mask (first index horizontal, second vertical)
    if ~any(m(:)), return; end
    C = zeros([size(m') 3]); for c = 1:3, C(:,:,c) = rgb(c); end
    h = image(C); set(h,'AlphaData',alpha*double(m'));
end
