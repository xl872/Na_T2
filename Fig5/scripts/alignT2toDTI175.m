%% =====================================================================
%  [175 version: same session as the sodium data -> useHeaderOnly = true, see settings]
%  Align the sodium T2* maps (Na23_scan142-160, 23Na scan 142 grid, 16^3 @ 0.40 mm)
%  to the proton DTI data (DWI_scan175, 1H scan 175 grid)
%  and compare ADC vs T2* (fast / slow / mono) voxel by voxel.
%
%  Run this from the DWI_scan175 folder, cell by cell (Ctrl+Enter).
%
%  REVISION 2026-09-21 (same changes as alignT2toDTI613New.m, but NO manual shift: same session, header-only alignment)
%    * T2estimate_scan142_160.m now fits the whole matrix (mask = allave > 12000), so the DWI ROI alone defines the
%      organoid area for ADC and for T2* (useSodiumROI = false) and the three maps of 6b share ONE mask
%    * fit bounds are read from T2estimate.m (readT2Bounds.m); fast / slow scatter drop only their own bound voxels
%    * the ROI analysis that follows this script is ROI_T2long_T2vsADC.m (ROIs = whole sodium voxels from the T2* long map)
%    * trimTopRows = 1: the top row of the DWI ROI (rim) is removed in every slice, for ADC and both T2* maps alike
%
%  HOW THE ALIGNMENT WORKS
%  -----------------------
%  Both AFNI grids came out of to3d with generic axis labels, so the +orig frames of
%  the two datasets are NOT a common physical frame (the sample was also moved between
%  the 23Na session, 2025-08-06, and the 1H session, 2025-09-11).  The alignment is
%  therefore done in two parts:
%    (a) ORIENTATION from the Bruker headers (pdata/1/visu_pars): voxel index -> scanner
%        (subject) mm for each dataset.  This fixes the axis permutations, the reversed
%        slice storage of the sodium scan (disk_reverse_slice_order = a mirror if ignored)
%        and the oblique slice orientation.  No sample-position information is trusted.
%    (b) RESIDUAL RIGID MOVE from image features: the sample tube (a cylinder -> its axis,
%        fitted from the disc centres of the coronal slices) and the organoid (ROI
%        centroid).  Tilt = rotate the sodium tube axis onto the proton tube axis,
%        Shift = make the organoid centroids coincide.  A roll about the tube axis cannot
%        be measured from a cylinder + a (near-axis) round organoid; the header roll is
%        kept (useRoll = false) unless the organoid sits clearly off-axis in both scans.
%  Aliasing / motion in 175: the tube is longer than the phase-encoded FOV along
%  the tube axis (z), so medium signal wraps in along z.  Only disc CENTRES per coronal
%  slice are used for the tube fit (the wrapped medium is uniform, it does not move the
%  centres) and nothing is taken from the tube ends.  Edge voxels of the proton ROI are
%  eroded before the ADC is box-averaged onto the sodium grid (erodeROI), and only sodium
%  voxels well covered by the proton ROI enter the statistics (minFrac).
%
%  OUTPUTS
%    DWI_scan175 (proton grid): T2_on175, T2s_on175, T2f_on175, ROI_T2_on175,
%                         allave_T2_on175 (+orig), T2_to_175_xform.mat,
%                         ADCvsT2_175_voxels.csv, QC_*.jpg, ADCvsT2_*.jpg
%    Na23_scan142-160 (sodium grid): ADC_175_on142grid, ROIfrac_175_on142grid (+orig)
%                         (next to the sodium allave so AFNI can overlay them)
%
%  KNOWN LIMITS
%    * roll about the tube axis is not measurable (see above)
%    * the per-voxel sodium biexponential fit pins most voxels at its bounds (T2f = 0.1 or 2 ms,
%      T2s = 5 or 20 ms, see T2estimate.m); exclBounds drops those from the fast/slow scatter
%    * the wrapped (aliased) medium signal adds a fast-decaying component to every 175 voxel,
%      including the organoid; that biases the per-voxel ADC upward mainly at the lowest b
%    * the same wrapped medium (laterally displaced because the tube is tilted) bends the proton
%      disc-centre track, so the proton tube-axis DIRECTION is uncertain by ~10 deg.  Because the
%      alignment is anchored on the organoid centroid, a 10 deg tilt error moves organoid voxels
%      by < 0.15 mm (half a sodium voxel); it does not matter for the ADC-vs-T2 comparison
%    * step 8 (optional) fits the sodium decay in ADC-defined zones (low/mid/high ADC) instead of per
%      voxel; that is the robust way to get fast/slow T2* vs ADC from these data.  It shows WHY the
%      per-voxel fits pin: the short component is ~2.5-3 ms (above the 2 ms fast bound of
%      T2estimate.m) and ~55 % of the signal is a long component
%      (absorbed by the constant term, so the "slow" component collapses onto its 5 ms bound).
%      Re-running T2estimate.m with e.g. fast [0.1 10] ms and slow [5 500] ms would give usable maps.
%  =====================================================================
cd(fullfile(fileparts(mfilename('fullpath')),'..','data','DWI_scan175'));   % release layout: run in the dataset folder
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','helper'));   % release layout: readT2Bounds.m
clc; clear; close all;

%% ---------------- 0. settings ----------------------------------------------------
DTI = pwd;                                          % .../DWI_scan175
T2D = fullfile(fileparts(pwd),'Na23_scan142-160');       % .../Na23_scan142-160
RAW = fullfile('..','bruker','20250801_201306_Mouse_23Naplus1_08012025_1_81');   % release layout
scanD = 175;   % proton DtiStandard (this folder)
scanT = 142;   % sodium scan that defines the Na23_scan142-160 grid (T2 series 142-160, 16^3 @ 0.40 mm)

useRoll       = false;     % also roll about the tube axis so the organoid's off-axis direction matches
shiftMode     = 'organoid';% 'organoid': organoid centroids coincide (default)
                           % 'tube'    : tube axes coincide, the organoid only fixes the position along the axis
                           % 'both'    : half way between the two
extraShift_mm = [0 0 0];   % manual tweak of the sodium data, mm along proton x(R-L) y(S-I) z(A-P); rerun 4-6
useHeaderOnly = true;      % scan 175 is from the SAME session as the sodium series 142-160 (2025-08-04/05) with the identical
                           % prescription (same VisuCoreOrientation, FOV 6.4 mm, position 0.1 mm apart): trust the headers, no
                           % feature-based tilt/shift.  The feature-based correction is still computed and printed as a check.
erodeROI      = 1;         % erode the proton ROI by n voxels before averaging ADC into sodium voxels
minFrac       = 0.5;       % sodium voxel must be covered >= this fraction by the (eroded) proton ROI
useSodiumROI  = false;     % false: the DWI ROI (roiD) alone defines the organoid area for ADC AND T2* (same area in both);
                           % true : additionally require the sodium ROI (Na23_scan142-160/ROI+orig), as before 2026-09-21
exclBounds    = true;      % fast / slow scatter: drop the sodium voxels whose T2f / T2s sits on its OWN fit bound
T2fBounds     = [];        % [] = read the per-voxel fit limits from T2estimate_scan142_160.m with readT2Bounds.m;
T2sBounds     = [];        %      or set by hand, e.g. [0 5] and [5 50] (ms)
trimTopRows   = 1;         % remove this many TOP rows (highest j = top of the coronal maps of 6b) of the DWI ROI in every slice:
                           % that row is the rim with partial-volume ADC, and the three maps then show the same, smaller area.
                           % The trimmed ROI is written as <ROI>_notop+orig and is the one used for the ADC box average (5b), the
                           % statistics (6a) and the maps (6b); ROI_T2long_T2vsADC.m reads it too.  0 = keep the whole DWI ROI
mapMargin     = 3;         % proton-grid maps (6b): show the organoid ROI box (ROI+orig, all slices) + this many proton
                           % voxels (0.16 mm) around it; [] = full FOV

fprintf('proton folder : %s\nsodium folder : %s\n',DTI,T2D);

%% ---------------- 1. load the AFNI datasets and have a look ----------------------
[anatD,InfoD] = BrikLoad(fullfile(DTI,'allave+orig.BRIK'));      anatD = double(anatD);
roiFileD = 'ROI_shrink';   % proton ROI used for the alignment (as in alignT2toDTI613New.m); falls back to ROI+orig when absent
if ~exist(fullfile(DTI,[roiFileD '+orig.BRIK']),'file'), fprintf('%s+orig not found in this folder -> using ROI+orig\n', roiFileD); roiFileD = 'ROI'; end
roiD  = BrikLoad(fullfile(DTI,[roiFileD '+orig.BRIK'])) > 0;
adcD  = double(BrikLoad(fullfile(DTI,'ADC+orig.BRIK')));          % ADC in 1e-6 mm^2/s (myADCProt.m)
[anatT,InfoT] = BrikLoad(fullfile(T2D,'allave+orig.BRIK'));      anatT = double(anatT);
roiT  = BrikLoad(fullfile(T2D,'ROI+orig.BRIK')) > 0;
T2m   = double(BrikLoad(fullfile(T2D,'T2+orig.BRIK')));            % mono-exponential T2* (ms)
T2s   = double(BrikLoad(fullfile(T2D,'T2s+orig.BRIK')));           % slow component (ms)
T2f   = double(BrikLoad(fullfile(T2D,'T2f+orig.BRIK')));           % fast component (ms)
fprintf('sodium maps   : T2 / T2s / T2f (%d voxels with a fit; %d of the %d sodium-ROI voxels)\n', nnz(T2m>0), nnz(T2m>0 & roiT), nnz(roiT));
if isempty(T2fBounds) || isempty(T2sBounds)
    [T2fBounds, T2sBounds, bnd] = readT2Bounds(fullfile('..','..','scripts','T2estimate_scan142_160.m'));
    fprintf('fit bounds read from T2estimate.m line %d: T2f [%g %g], T2s [%g %g] ms\n', bnd.line, T2fBounds, T2sBounds);
end
ND = size(anatD); NT = size(anatT);
roiM = roiD;                                   % DWI ROI used for maps / ADC averaging / statistics (roiD itself: alignment QC only)
for k = 1:ND(3)*(trimTopRows > 0)
    jj = find(any(roiD(:,:,k),1));             % rows (2nd index j, vertical in the maps) that hold ROI voxels in this slice
    if ~isempty(jj), roiM(:, jj(max(1,end-trimTopRows+1):end), k) = false; end
end
if trimTopRows > 0
    fprintf('DWI ROI: top %d row(s) removed in every slice: %d -> %d voxels\n', trimTopRows, nnz(roiD), nnz(roiM));
    writeOut(double(roiM), InfoD, fullfile(DTI,[roiFileD '_notop']));
end
resD = abs(InfoD.DELTA(1)); resT = abs(InfoT.DELTA(1));
fprintf('proton  %dx%dx%d @ %.2f mm, ROI %d voxels (%.2f mm^3)\n',ND,resD,nnz(roiD),nnz(roiD)*resD^3);
fprintf('sodium  %dx%dx%d @ %.2f mm, ROI %d voxels (%.2f mm^3)\n',NT,resT,nnz(roiT),nnz(roiT)*resT^3);

% organoid slice (3rd index = A-P = tube axis in both datasets, coronal = round view)
[iD,jD,kD] = ind2sub(ND,find(roiD)); cD_ijk = [mean(iD) mean(jD) mean(kD)];   % 1-based
[iT,jT,kT] = ind2sub(NT,find(roiT)); cT_ijk = [mean(iT) mean(jT) mean(kT)];
figure('Name','1 native data','Position',[50 50 1100 420]);
subplot(1,2,1); showSlice(anatD(:,:,round(cD_ijk(3))), roiD(:,:,round(cD_ijk(3))), 'r');
title(sprintf('proton allave, coronal k=%d (round view)  red = proton ROI',round(cD_ijk(3))));
subplot(1,2,2); showSlice(anatT(:,:,round(cT_ijk(3))), roiT(:,:,round(cT_ijk(3))), 'r');
title(sprintf('sodium allave, coronal k=%d (round view)  red = sodium ROI',round(cT_ijk(3))));

%% ---------------- 2. true geometry from the Bruker headers -----------------------
%  A * [i;j;k;1] (0-based AFNI indices) = [x;y;z;1] subject mm (+x=L, +y=P, +z=S).
%  proton 175   : AFNI (i,j,k) = 2dseq (d1, d2, 80-d3_disk)  (proc175_to3d.sh: true labels R-L/S-I/P-A + 3daxialize RSA; slices stored in reverse order)
%  sodium 142    : AFNI (i,j,k) = 2dseq (d1, d2, d3)      (generic to3d labels, no permutation)
PD = eye(4);                                            % DWI_scan175: generic to3d labels, AFNI (i,j,k) = 2dseq (d1,d2,d3); verified against the raw 2dseq (corr 1.000)
PT = eye(4);                                            % ijk -> d1 d2 d3 for the sodium
[A_D,hdrD] = brukerIJK2subject(fullfile(RAW,num2str(scanD),'pdata','1','visu_pars'), PD);
[A_T,hdrT] = brukerIJK2subject(fullfile(RAW,num2str(scanT),'pdata','1','visu_pars'), PT);
describeFrame('proton 175', A_D, ND, hdrD);
describeFrame('sodium 142', A_T, NT, hdrT);
% sanity: voxel size recovered from the affine must equal the AFNI voxel size
assert(abs(norm(A_D(1:3,1))-resD)<1e-6 && abs(norm(A_T(1:3,1))-resT)<1e-6,'header/AFNI voxel size mismatch');

%% ---------------- 3. features: tube axis and organoid, in subject mm --------------
% tube lumen (medium, bright in both) -> disc centre in every coronal slice -> 3D line fit
[lumD, ctrD, okD] = tubeDiscCentres(anatD);     % centres are 0-based (i,j,k) per slice
[lumT, ctrT, okT] = tubeDiscCentres(anatT);
XD = toSubject(A_D, ctrD(okD,:));  XT = toSubject(A_T, ctrT(okT,:));
[tcD,dD] = fitLine(XD);  [tcT,dT] = fitLine(XT);
if dD(2) < 0, dD = -dD; end                       % point the axes the same way (+y = P)
if dot(dT,dD) < 0, dT = -dT; end
% organoid centroid (ROI) in subject mm
cD = toSubject(A_D, cD_ijk-1)';  cT = toSubject(A_T, cT_ijk-1)';     % columns, subject mm
pD = radialOffset(cD,tcD,dD);    pT = radialOffset(cT,tcT,dT);
ang = acosd(min(1,dot(dD,dT)));
fprintf('\nTube axis (subject frame)  proton [%6.3f %6.3f %6.3f]  sodium [%6.3f %6.3f %6.3f]  -> %.1f deg apart\n',dD,dT,ang);
fprintf('Tube fit residual (rms, mm)  proton %.3f (%d/%d slices)  sodium %.3f (%d/%d slices)\n', ...
        lineResidual(XD,tcD,dD),nnz(okD),ND(3),lineResidual(XT,tcT,dT),nnz(okT),NT(3));
fprintf('Organoid radial offset from tube axis  proton %.2f mm  sodium %.2f mm  (should be similar; >~0.5 mm needed for useRoll)\n',norm(pD),norm(pT));
fprintf('Tube diameter estimate  proton %.2f mm  sodium %.2f mm\n', 2*sqrt(mean(sum(lumD,[1 2]),'all')/pi)*resD, 2*sqrt(mean(sum(lumT,[1 2]),'all')/pi)*resT);
% how stable is the proton axis direction?  fit the two halves of the FOV separately (aliased medium
% and B1 shading bend the thresholded disc centres a little)
hD = ctrD(okD,3) < median(ctrD(okD,3)); [~,d1] = fitLine(XD(hD,:)); [~,d2] = fitLine(XD(~hD,:));
fprintf('Proton axis from the two halves of the FOV differs by %.1f deg (uncertainty of the tilt)\n', acosd(min(1,abs(dot(d1,d2)))));

figure('Name','3 tube fit','Position',[50 500 1100 420]);
subplot(1,2,1); plot(ctrD(:,3),ctrD(:,1),'r.-',ctrD(:,3),ctrD(:,2),'b.-'); hold on
plot(ctrD(okD,3),ctrD(okD,1),'ro',ctrD(okD,3),ctrD(okD,2),'bo'); xlabel('k (A-P slice)'); ylabel('disc centre (voxel)');
legend('i (R-L)','j (S-I)','used','used','Location','best'); title('proton: tube disc centre vs slice (tilt of the tube in the FOV)'); grid on
subplot(1,2,2); plot(ctrT(:,3),ctrT(:,1),'r.-',ctrT(:,3),ctrT(:,2),'b.-'); hold on
plot(ctrT(okT,3),ctrT(okT,1),'ro',ctrT(okT,3),ctrT(okT,2),'bo'); xlabel('k (A-P slice)'); ylabel('disc centre (voxel)');
legend('i (R-L)','j (S-I)','used','used','Location','best'); title('sodium: tube disc centre vs slice'); grid on
print(gcf,fullfile(DTI,'QC_tubefit_175.jpg'),'-djpeg','-r150');

%% ---------------- 4. rigid correction  Q : sodium subject -> proton subject --------
Rtilt = rotAtoB(dT,dD);                              % tube axes parallel
pT2   = Rtilt*pT;
phi   = atan2(dot(cross(pT2,pD),dD), dot(pT2,pD));   % roll that would align the organoid's off-axis direction
if useRoll, RQ = rotAxis(dD,phi)*Rtilt; else, RQ = Rtilt; end
shiftS = A_D(1:3,1:3)*(extraShift_mm(:)/resD);      % manual tweak, proton axes mm -> subject mm
t0  = cD - RQ*cT;                                    % translation that makes the organoid centroids coincide
rAx = radialOffset(RQ*tcT + t0, tcD, dD);           % ... and where the sodium tube axis then lands vs the proton axis
switch shiftMode
    case 'organoid', t = t0;
    case 'tube',     t = t0 - rAx;
    case 'both',     t = t0 - rAx/2;
    otherwise, error('shiftMode must be organoid, tube or both');
end
if useHeaderOnly
    fprintf('\nuseHeaderOnly: feature-based correction NOT applied (it would have been tilt %.1f deg, shift [%.2f %.2f %.2f] mm)\n', ang, t);
    RQ = eye(3); t = zeros(3,1);
end
Q = eye(4); Q(1:3,1:3) = RQ; Q(1:3,4) = t + shiftS;
M_T2ijk_to_Dijk = A_D \ (Q*A_T);                     % sodium ijk -> proton ijk (0-based)
M_Dijk_to_T2ijk = inv(M_T2ijk_to_Dijk);
fprintf('\nRigid correction: tilt %.1f deg, roll %s (%.1f deg would align the organoid azimuth)\n', ang, string(useRoll), rad2deg(phi));
fprintf('Organoid position by header alone differs by %.2f mm between the two sessions (= how far the sample was moved)\n', norm(cD-cT));
fprintf('Sodium -> proton translation in subject frame: [%.2f %.2f %.2f] mm  (shiftMode = %s)\n',Q(1:3,4),shiftMode);
fprintf('Consistency of the two landmarks: with organoid centroids matched the tube axes are %.2f mm apart\n', norm(rAx));
tcT_Q = RQ*tcT + Q(1:3,4);
fprintf('After this shift: tube axes %.2f mm apart (%.2f deg), organoid centroids %.2f mm apart\n', ...
        norm(radialOffset(tcT_Q,tcD,dD)), acosd(min(1,abs(dot(RQ*dT,dD)))), norm(RQ*cT+Q(1:3,4)-cD));
q = M_T2ijk_to_Dijk*[cT_ijk(:)-1;1];
fprintf('Check: sodium organoid centroid -> proton ijk = [%.1f %.1f %.1f] (proton ROI centroid [%.1f %.1f %.1f], 1-based)\n', q(1:3)'+1, cD_ijk);
save(fullfile(DTI,'T2_to_175_xform.mat'),'Q','A_D','A_T','M_T2ijk_to_Dijk','M_Dijk_to_T2ijk','dD','dT','cD','cT','tcD','tcT','useRoll','shiftMode','extraShift_mm','useHeaderOnly');

%% ---------------- 5. resample both ways and write AFNI datasets --------------------
% 5a. sodium -> proton grid (nearest neighbour; one sodium voxel covers 3x3x3 proton voxels)
[I,J,K] = ndgrid(0:ND(1)-1,0:ND(2)-1,0:ND(3)-1);
tq  = M_Dijk_to_T2ijk*[I(:)';J(:)';K(:)';ones(1,numel(I))];   % continuous sodium ijk (0-based)
tr  = round(tq(1:3,:));
in  = all(tr>=0,1) & tr(1,:)<=NT(1)-1 & tr(2,:)<=NT(2)-1 & tr(3,:)<=NT(3)-1;
lin = sub2ind(NT,tr(1,in)+1,tr(2,in)+1,tr(3,in)+1);
pull = @(v) reshape(accum0(v,lin,in,numel(I)),ND);
T2m_on = pull(T2m); T2s_on = pull(T2s); T2f_on = pull(T2f); roiT_on = pull(double(roiT))>0.5; lumT_on = pull(double(lumT))>0.5;
anatT_on = zeros(ND); anatT_on(in) = interp3(anatT, tq(2,in)+1, tq(1,in)+1, tq(3,in)+1, 'linear', 0);  % (x<->dim2, y<->dim1 in interp3)
writeOut(T2m_on ,InfoD,fullfile(DTI,'T2_on175'));
writeOut(T2s_on ,InfoD,fullfile(DTI,'T2s_on175'));
writeOut(T2f_on ,InfoD,fullfile(DTI,'T2f_on175'));
writeOut(double(roiT_on),InfoD,fullfile(DTI,'ROI_T2_on175'));
writeOut(anatT_on,InfoD,fullfile(DTI,'allave_T2_on175'));

% 5b. proton ADC -> sodium grid (box average of the proton ROI voxels that fall in each sodium voxel)
roiDe = roiM; for e = 1:erodeROI, roiDe = imerode(roiDe, strel('sphere',1)); end
rv = double(roiDe(in)); rv = rv(:);  av = adcD(in); av = av(:);
cntAll = accumarray(lin(:),1,[prod(NT) 1]);                          % proton voxels per sodium voxel (~27)
cntRoi = accumarray(lin(:),rv,[prod(NT) 1]);                         % of which inside the eroded proton ROI
sumAdc = accumarray(lin(:),rv.*av,[prod(NT) 1]);
ADC_onT  = reshape(sumAdc./max(cntRoi,1),NT);  ADC_onT(cntRoi==0) = 0;
frac_onT = reshape(cntRoi./max(cntAll,1),NT);
writeOut(ADC_onT ,InfoT,fullfile(T2D,'ADC_175_on142grid'));       % sodium grid -> next to the sodium allave
writeOut(frac_onT,InfoT,fullfile(T2D,'ROIfrac_175_on142grid'));

% 5c. QC: overlap of the two organoid ROIs and of the tube edges
dice = 2*nnz(roiD & roiT_on)/(nnz(roiD)+nnz(roiT_on));
[i2,j2,k2] = ind2sub(ND,find(roiT_on)); cT_onD = [mean(i2) mean(j2) mean(k2)];
fprintf('\nQC on proton grid: ROI Dice = %.2f, centroid distance = %.2f mm, sodium ROI mapped = %d voxels\n', ...
        dice, norm((cT_onD-cD_ijk)*resD), nnz(roiT_on));
fprintf('QC on sodium grid: %d sodium voxels touched by the (eroded) proton ROI, %d with coverage >= %.2f, %d of them inside the sodium ROI\n', ...
        nnz(frac_onT>0), nnz(frac_onT>=minFrac), minFrac, nnz(frac_onT>=minFrac & roiT));
kk = round(cD_ijk(3)); jj = round(cD_ijk(2)); ii = round(cD_ijk(1)); kT = round(cT_ijk(3));
% tube-vs-tube check on the proton grid: disc centres of the proton lumen and of the mapped sodium lumen, slice by slice
offs = nan(ND(3),1); dP = nan(ND(3),1); dS = nan(ND(3),1);
for k = 1:ND(3)
    a = lumD(:,:,k); b = lumT_on(:,:,k);
    if nnz(a) > 100 && nnz(b) > 100
        [ia,ja] = find(a); [ib,jb] = find(b);
        offs(k) = norm([mean(ia)-mean(ib), mean(ja)-mean(jb)])*resD; dP(k) = 2*sqrt(nnz(a)/pi)*resD; dS(k) = 2*sqrt(nnz(b)/pi)*resD;
    end
end
fprintf('Tube overlap on proton grid (%d slices): disc-centre offset %.2f +/- %.2f mm, diameters proton %.2f mm / sodium %.2f mm\n', ...
        nnz(~isnan(offs)), mean(offs,'omitnan'), std(offs,'omitnan'), mean(dP,'omitnan'), mean(dS,'omitnan'));
figure('Name','5 alignment QC','Position',[50 50 1500 420]);
subplot(1,4,1); showSlice(anatD(:,:,kk), roiD(:,:,kk),'r', roiT_on(:,:,kk),'c', lumT_on(:,:,kk),'y');
title(sprintf('proton coronal k=%d: red proton ROI, cyan sodium ROI, yellow sodium tube edge',kk),'FontSize',8);
subplot(1,4,2); showSlice(squeeze(anatD(:,jj,:)), squeeze(roiD(:,jj,:)),'r', squeeze(roiT_on(:,jj,:)),'c', squeeze(lumT_on(:,jj,:)),'y');
title(sprintf('proton axial j=%d (x horiz, z tube axis vert)',jj),'FontSize',8);
subplot(1,4,3); showSlice(squeeze(anatD(ii,:,:)), squeeze(roiD(ii,:,:)),'r', squeeze(roiT_on(ii,:,:)),'c', squeeze(lumT_on(ii,:,:)),'y');
title(sprintf('proton sagittal i=%d (y horiz, z vert)',ii),'FontSize',8);
subplot(1,4,4); showSlice(anatT(:,:,kT), roiT(:,:,kT),'r', frac_onT(:,:,kT)>=minFrac,'c');
title(sprintf('sodium coronal k=%d: red sodium ROI, cyan proton ROI (cov>=%.1f)',kT,minFrac),'FontSize',8);
print(gcf,fullfile(DTI,'QC_alignment_175.jpg'),'-djpeg','-r150');

%% ---------------- 6. voxel-wise ADC vs T2* ------------------------------------------
% 6a. on the sodium grid (independent voxels; each has its own T2* fit and the mean ADC of the
%     proton ROI voxels inside it)  -> use this one for statistics
sel = frac_onT>=minFrac & ADC_onT>0 & T2m>0;          % the DWI ROI defines the area; the whole-matrix T2* maps have a fit everywhere
if useSodiumROI, sel = sel & roiT; end
tol = 1e-3;
onBf = abs(T2f-T2fBounds(1))<tol | abs(T2f-T2fBounds(2))<tol;      % T2f on its own bound
onBs = abs(T2s-T2sBounds(1))<tol | abs(T2s-T2sBounds(2))<tol;      % T2s on its own bound
atBound = onBf | onBs;  ownB = {onBf,onBs,false(NT)};
fprintf('\nSodium-grid comparison: %d voxels (%d of them inside the sodium ROI); on a fit bound: T2f %d, T2s %d (own-bound voxels dropped per panel: %s)\n', ...
        nnz(sel), nnz(sel&roiT), nnz(sel&onBf), nnz(sel&onBs), string(exclBounds));
figure('Name','6 ADC vs T2 (sodium grid)','Position',[50 50 1200 380]);
names = {'T2* fast (ms)','T2* slow (ms)','T2* mono (ms)'}; maps = {T2f,T2s,T2m};
for p = 1:3
    selp = sel; if exclBounds, selp = sel & ~ownB{p}; end
    adc = ADC_onT(selp); v = maps{p}(selp);
    subplot(1,3,p); scatter(adc,v,30,'filled'); hold on
    if nnz(selp) >= 3
        lsline; [rs,ps] = corr(adc,v,'type','Spearman'); [rp,pp] = corr(adc,v,'type','Pearson');
        title(sprintf('n=%d  Spearman \\rho=%.2f (p=%.3f)  Pearson r=%.2f (p=%.3f)',nnz(selp),rs,ps,rp,pp),'FontSize',8);
        fprintf('  ADC vs %-14s  n=%3d  Spearman rho=%+.2f p=%.3f   Pearson r=%+.2f p=%.3f\n',names{p},nnz(selp),rs,ps,rp,pp);
    else
        title(sprintf('n=%d (too few voxels)',nnz(selp)),'FontSize',8);
    end
    xlabel('ADC (um^2/s)'); ylabel(names{p});
end
sgtitle(sprintf('ADC vs sodium T2*, per sodium voxel (ROI coverage >= %.1f, proton ROI eroded %d vox)',minFrac,erodeROI));
print(gcf,fullfile(DTI,'ADCvsT2_175_sodiumgrid.jpg'),'-djpeg','-r150');
[si,sj,sk] = ind2sub(NT,find(sel));
writetable(table(si,sj,sk,ADC_onT(sel),T2f(sel),T2s(sel),T2m(sel),frac_onT(sel),atBound(sel), ...
           'VariableNames',{'i','j','k','ADC_1e-6mm2s','T2fast_ms','T2slow_ms','T2mono_ms','ROIcoverage','fitOnBound'}), ...
           fullfile(DTI,'ADCvsT2_175_voxels.csv'));
fprintf('  wrote ADCvsT2_175_voxels.csv\n');

%% 6b. on the proton grid (0.16 mm; T2* is nearest-neighbour, i.e. constant over ~2.5^3 blocks) -> maps only
%  ONE mask for the three maps = the DWI ROI (roiD): ADC and T2* cover exactly the same area (count printed).
selD = roiM & T2s_on>0 & T2f_on>0 & adcD>0;
kk=21;   % coronal slice shown (the slice of page 3 of Organoid_supp.pptx; myADCProt.m shows 20); alignT2toDTI613New.m uses kk=14
fprintf('\nproton-grid maps, slice k=%d: DWI ROI %d voxels, with ADC and T2* %d (all slices: %d of %d)\n', kk, nnz(roiM(:,:,kk)), nnz(selD(:,:,kk)), nnz(selD), nnz(roiM));
roiBox = BrikLoad(fullfile(DTI,'ROI+orig.BRIK')) > 0;     % main organoid ROI: defines the field of view of the maps (mapMargin)
figure('Name','6 maps on proton grid','Position',[50 500 1500 400]);
subplot(1,3,1); showMap(anatD(:,:,kk), adcD(:,:,kk), selD(:,:,kk), 'ADC (um^2/s)',[100 300]); zoomToROI(roiBox,mapMargin); title(sprintf('ADC, coronal k=%d',kk));
subplot(1,3,2); showMap(anatD(:,:,kk), T2s_on(:,:,kk), selD(:,:,kk), 'T2* long (ms)',[0 15]); zoomToROI(roiBox,mapMargin); title('sodium T2* slow on proton grid');
subplot(1,3,3); showMap(anatD(:,:,kk), T2f_on(:,:,kk), selD(:,:,kk), 'T2* short (ms)',[0 1.5]); zoomToROI(roiBox,mapMargin); title('sodium T2* fast on proton grid');
set(findall(gcf,'Type','axes'),'Toolbar',[]);              % no axes-toolbar icons in a batch-mode export
print(gcf,fullfile(DTI,['ADCvsT2_175_protongrid_slice',num2str(kk),'.jpg']),'-djpeg','-r150');

%% ---------------- 7. optional: did the sample move during the 175 scan? -------------
% organoid centre in each b-value volume (allave_t).  A drift of >1 voxel between frames
% means the per-voxel ADC fit mixes neighbouring voxels -> trust the eroded/box-averaged values more.
dat = double(BrikLoad(fullfile(DTI,'allave_t+orig.BRIK')));
dro = imdilate(roiD, strel('sphere',3));
fprintf('\nOrganoid centre per b-value volume (voxels, 1-based):\n');
for t = 1:size(dat,4)
    v = dat(:,:,:,t); thr = prctile(v(dro),60); m = v>thr & dro;
    [a,b,c] = ind2sub(ND,find(m)); fprintf('  frame %d : i=%.2f j=%.2f k=%.2f\n',t,mean(a),mean(b),mean(c));
end

%% ---------------- 8. optional: sodium decay in ADC-defined zones ---------------------
%  Per-voxel biexponential fits of the sodium data mostly end on their bounds (too little SNR per
%  voxel and bounds that do not match this sample).  Here the organoid (sodium voxels covered
%  by the proton ROI) is split into three zones by the proton ADC mapped onto the sodium grid, the raw
%  sodium signal is averaged per zone over all TEs (same runs / smoothing as T2estimate.m) and fitted
%  once per zone with WIDER bounds: fast 0.1-10 ms, slow 5-500 ms.
runs  = 142:160;                                              % Na23_scan142-160 T2 series (as in T2estimate_scan142_160.m)
dname = '20250801_201306_Mouse_23Naplus1_08012025_1_81';
TE = zeros(size(runs)); S = zeros([NT numel(runs)]);
for r = 1:numel(runs)
    v = double(BrikLoad(fullfile(T2D,sprintf('%s.%d.an200epi+orig.BRIK',dname,runs(r)))));
    if ndims(v) == 4, v = mean(v,4); end
    S(:,:,:,r) = imgaussfilt3(v,0.5);
    TE(r) = brukerParam(fullfile(RAW,num2str(runs(r)),'method'),'PVM_EchoTime');   % ms
end
sel8 = frac_onT>=minFrac & ADC_onT>0;  if useSodiumROI, sel8 = sel8 & roiT; end
qs = quantile(ADC_onT(sel8),[1/3 2/3]);
zone = zeros(NT); zone(sel8 & ADC_onT<=qs(1)) = 1; zone(sel8 & ADC_onT>qs(1) & ADC_onT<=qs(2)) = 2; zone(sel8 & ADC_onT>qs(2)) = 3;
zname = {'low ADC','mid ADC','high ADC'}; col = {'b','g','r'};
fprintf('\nSodium decay per ADC zone (%d voxels, ADC tertile cuts %.0f / %.0f um^2/s):\n',nnz(sel8),qs);
figure('Name','8 decay per ADC zone','Position',[50 50 1000 420]);
Sv = reshape(S,[],numel(runs));
for z = 1:3
    m = zone==z; y = mean(Sv(m(:),:),1)';  x = TE(:);
    if nnz(m) < 2 || any(~isfinite(y)), fprintf('  %-8s n=%2d  -> too few voxels for a zone fit, skipped\n', zname{z}, nnz(m)); continue; end
    % biexponential, no constant: a*((1-d) exp(-x/b) + d exp(-x/e)), wide bounds
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0.1 0 5],'Upper',[Inf 10 1 500],'StartPoint',[y(1) 2 0.5 50]);
    fb = fit(x,y,fittype('a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo));
    % short component + constant (= T2estimate.m mono model): b = short T2*, c/(a+c) = long fraction
    fo1 = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],'Upper',[Inf 50 Inf],'StartPoint',[y(1) 3 y(end)]);
    fm = fit(x,y,fittype('c+a*exp(-x/b)','options',fo1));
    fprintf('  %-8s n=%2d  mean ADC %5.0f | biexp: T2*fast %.2f ms  T2*slow %.0f ms  long fraction %.2f | short+const: T2*short %.2f ms  long fraction %.2f\n', ...
            zname{z}, nnz(m), mean(ADC_onT(m)), fb.b, fb.e, fb.d, fm.b, fm.c/(fm.a+fm.c));
    subplot(1,2,1); semilogy(x,y/y(1),[col{z} 'o']); hold on; xx = linspace(min(x),max(x),200)'; semilogy(xx,fb(xx)/y(1),[col{z} '-']);
    subplot(1,2,2); bar(z-0.2,fb.b,0.35,'FaceColor',col{z}); hold on; bar(z+0.2,fm.b,0.35,'FaceColor',col{z},'FaceAlpha',0.4);
end
subplot(1,2,1); xlabel('TE (ms)'); ylabel('sodium signal (normalised to first TE)'); legend(reshape([zname; repmat({''},1,3)],1,[]),'Location','best'); title('sodium decay per ADC zone (biexp fit, wide bounds)'); grid on
subplot(1,2,2); set(gca,'XTick',1:3,'XTickLabel',zname); ylabel('short T2* (ms)'); title('short T2* per zone (solid: biexp fast, pale: short+const)');
print(gcf,fullfile(DTI,'ADCzones_sodiumdecay_175.jpg'),'-djpeg','-r150');

%% ================= local functions ===================================================
function [A,hdr] = brukerIJK2subject(visuFile, Pmap)
% affine 0-based AFNI ijk -> subject mm from visu_pars; Pmap maps AFNI ijk -> 2dseq (d1,d2,d3)
    hdr.size   = brukerParam(visuFile,'VisuCoreSize');
    hdr.extent = brukerParam(visuFile,'VisuCoreExtent');
    hdr.orient = reshape(brukerParam(visuFile,'VisuCoreOrientation'),3,3)';   % rows = direction of d1,d2,d3
    hdr.pos    = brukerParam(visuFile,'VisuCorePosition');   hdr.pos = hdr.pos(1:3);
    hdr.order  = brukerParam(visuFile,'VisuCoreDiskSliceOrder');
    res = hdr.extent(:)'./hdr.size(:)';  hdr.res = res;
    B = eye(4);  B(1:3,1:3) = (hdr.orient'.*res);  B(1:3,4) = hdr.pos(:);     % geometric d -> mm
    F = eye(4);
    if contains(hdr.order,'reverse')       % slice k on disk is geometric slice N3-1-k
        F(3,3) = -1; F(3,4) = hdr.size(3)-1;
    end
    A = B*F*Pmap;
end

function v = brukerParam(file, name)
    txt = fileread(file);
    tok = regexp(txt,['##\$' name '=([^\n]*)\n([^#]*)'],'tokens','once');
    if isempty(tok), error('parameter %s not found in %s',name,file); end
    if startsWith(strtrim(tok{1}),'(')                      % array: values on following lines
        v = sscanf(regexprep(tok{2},'<|>',' '),'%f');
        if isempty(v), v = strtrim(regexprep(tok{2},'<|>','')); end
        if isnumeric(v), v = v(:)'; end
    else
        v = strtrim(tok{1});  num = str2double(v); if ~isnan(num), v = num; end
    end
end

function describeFrame(name, A, N, hdr)
    ax = {'i (x)','j (y)','k (z)'}; lab = {'L','P','S'}; neg = {'R','A','I'};
    fprintf('\n%s: %dx%dx%d, 2dseq res [%.2f %.2f %.2f] mm, disk order: %s\n',name,N,hdr.res,hdr.order);
    for a = 1:3
        d = A(1:3,a)/norm(A(1:3,a)); [~,m] = max(abs(d)); s = lab{m}; if d(m)<0, s = neg{m}; end
        fprintf('   AFNI %s runs towards subject %s   dir = [%6.3f %6.3f %6.3f]\n',ax{a},s,d);
    end
    c = A*[(N(:)-1)/2;1]; fprintf('   FOV centre in subject frame = [%.2f %.2f %.2f] mm (header only, sample was moved)\n',c(1:3));
end

function X = toSubject(A, ijk0)          % rows of 0-based ijk -> rows of subject mm
    X = (A*[ijk0'; ones(1,size(ijk0,1))])'; X = X(:,1:3);
end

function [lum, ctr, ok] = tubeDiscCentres(vol)
% threshold (Otsu, organoid brightness clipped) -> per coronal slice: fill, largest blob, centre (0-based)
    hi  = prctile(vol(:),97); a = min(vol,hi)/hi; thr = graythresh(a)*hi;
    n3 = size(vol,3); lum = false(size(vol)); ctr = nan(n3,3); area = zeros(n3,1);
    for k = 1:n3
        s = imfill(vol(:,:,k) > thr,'holes'); s = bwareafilt(s,1);
        lum(:,:,k) = s; area(k) = nnz(s);
        [ii,jj] = find(s); ctr(k,:) = [mean(ii)-1, mean(jj)-1, k-1];
    end
    ok = area > 0.8*median(area(area>0));      % drop slices where the disc is clipped / broken
    % also drop slices whose disc touches the image border (clipped disc -> biased centre)
    for k = 1:n3
        s = lum(:,:,k); if any(s(1,:)) || any(s(end,:)) || any(s(:,1)) || any(s(:,end)), ok(k) = false; end
    end
    if nnz(ok) < 3, ok = area > 0.8*median(area(area>0)); end
end

function [c,d] = fitLine(X)              % least-squares 3D line through points (rows)
    c = mean(X,1); [~,~,V] = svd(X-c,0); d = V(:,1); c = c(:);
end
function r = lineResidual(X,c,d)
    P = (X-c'); P = P - (P*d)*d'; r = sqrt(mean(sum(P.^2,2)));
end
function p = radialOffset(x,c,d)          % component of (x-c) perpendicular to d
    x = x(:); c = c(:); d = d(:); v = x-c; p = v - dot(v,d)*d;
end
function R = rotAtoB(a,b)                 % rotation taking unit vector a onto b
    a = a(:)/norm(a); b = b(:)/norm(b); v = cross(a,b); s = norm(v); c = dot(a,b);
    if s < 1e-12, R = eye(3); return; end
    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]; R = eye(3) + vx + vx*vx*((1-c)/s^2);
end
function R = rotAxis(u,ang)               % rotation by ang about unit axis u
    u = u(:)/norm(u); K = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]; R = eye(3) + sin(ang)*K + (1-cos(ang))*(K*K);
end
function out = accum0(v,lin,in,n)         % pull v(lin) into a vector of length n at positions in
    out = zeros(n,1); out(in) = v(lin);
end
function writeOut(M, Info, prefixPath)   % write a float dataset on the grid of Info
    [p,f] = fileparts(prefixPath); old = dir(fullfile(p,[f '+orig.*'])); for q = 1:numel(old), delete(fullfile(p,old(q).name)); end
    Info.BRICK_TYPES = 3; Info.DATASET_RANK(2) = 1; Info.BRICK_FLOAT_FACS = 0;
    Info.BRICK_STATS = [min(M(:)) max(M(:))];
    for fld = {'TAXIS_NUMS','TAXIS_FLOATS','TAXIS_OFFSETS','BRICK_STATAUX','BRICK_KEYWORDS'}
        if isfield(Info,fld{1}), Info = rmfield(Info,fld{1}); end
    end
    if isfield(Info,'BRICK_LABS'), Info.BRICK_LABS = f; end
    Opt.Prefix = fullfile(p,f); Opt.OverWrite = 'y'; Opt.verbose = 0; Opt.Scale = 0; Opt.View = '+orig';
    [err,ErrMessage] = WriteBrik(double(M),Info,Opt); if err, error(ErrMessage); end
    fprintf('  wrote %s+orig\n',f);
end
function showSlice(img, varargin)         % image slice (i horizontal? no: i vertical->use transpose) with contours
    hi = prctile(img(:),96); if ~(hi>0), hi = 1; end       % 96th pct: tube medium visible, organoid saturates
    imagesc(img'); axis image xy off; colormap(gca,'gray'); caxis([0 hi]); hold on
    for c = 1:2:numel(varargin)
        m = varargin{c}; if any(m(:)), contour(double(m'),[0.5 0.5],varargin{c+1},'LineWidth',1); end
    end
end
function showMap(anat, map, mask, lab,cclim)
    B = repmat(anat'/prctile(anat(:),99.5),1,1,3); B(B>1) = 1; image(B); axis image xy off; hold on
    h = imagesc(map'); set(h,'AlphaData',double(mask')); colormap(gca,'jet');
    vals = map(mask); if any(vals), caxis([prctile(vals,2) prctile(vals,98)]); end
    hb = colorbar; ylabel(hb,lab);
    clim([cclim(1) cclim(2)])
end
function zoomToROI(roi, margin)           % crop the current map to the ROI box (over all slices) + margin whole voxels; [] = full FOV
    if isempty(margin) || ~any(roi(:)), return; end
    [bi,bj] = ind2sub([size(roi,1) size(roi,2)], find(any(roi,3)));
    xlim([max(1,min(bi)-margin)-0.5  min(size(roi,1),max(bi)+margin)+0.5]);
    ylim([max(1,min(bj)-margin)-0.5  min(size(roi,2),max(bj)+margin)+0.5]);
end
