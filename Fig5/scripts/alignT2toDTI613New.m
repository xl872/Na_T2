%% =====================================================================
%  Align the sodium T2* maps (Na23_scan402-418, 23Na scan 402 grid)
%  to the proton DTI data (DWI_scan613, 1H scan 613 grid, 80x64x64 @ 0.10 mm)
%  and compare ADC vs T2* (fast / slow / mono) voxel by voxel.
%
%  Run this from the DWI_scan613 folder, cell by cell (Ctrl+Enter).
%
%  REVISION 2026-09-21
%    * extraShift_mm = [-0.1 0 0]: the sodium maps are moved 0.1 mm to the left of the coronal maps so that the central
%      highest T2* long voxels sit on the high-ADC island (step 4b prints the evidence for every shift)
%    * T2set = 'all': the whole-matrix sodium maps T2all / T2alls / T2allf are used, so the DWI ROI alone defines the
%      organoid area for ADC and for T2* (useSodiumROI = false) and the three maps of 6b share ONE mask
%    * fit bounds are read from T2estimate.m (readT2Bounds.m); fast / slow scatter drop only their own bound voxels
%    * the ROI analysis that follows this script is ROI_T2long_T2vsADC.m (ROIs = whole sodium voxels from the T2* long map)
%
%  HOW THE ALIGNMENT WORKS
%  -----------------------
%  Both AFNI grids came out of to3d with generic axis labels, so the +orig frames of
%  the two datasets are NOT a common physical frame (the sample was also moved between
%  the 23Na session, 2025-08-22, and the 1H session, 2025-09-08).  The alignment is
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
%  Aliasing / motion in 613: the tube is longer than the 6.4 mm phase-encode FOV along
%  the tube axis (z), so medium signal wraps in along z.  Only disc CENTRES per coronal
%  slice are used for the tube fit (the wrapped medium is uniform, it does not move the
%  centres) and nothing is taken from the tube ends.  Edge voxels of the proton ROI are
%  eroded before the ADC is box-averaged onto the sodium grid (erodeROI), and only sodium
%  voxels well covered by the proton ROI enter the statistics (minFrac).
%
%  OUTPUTS
%    DWI_scan613 (proton grid): T2_on613New, T2s_on613New, T2f_on613New, ROI_T2_on613New,
%                         allave_T2_on613New (+orig), T2_to_613New_xform.mat,
%                         ADCvsT2_613New_voxels.csv, QC_*.jpg, ADCvsT2_*.jpg
%    Na23_scan402-418 (sodium grid): ADC_613New_on402grid, ROIfrac_613New_on402grid (+orig)
%                         (next to the sodium allave so AFNI can overlay them)
%
%  KNOWN LIMITS
%    * roll about the tube axis is not measurable (see above)
%    * the per-voxel sodium biexponential fit pins most voxels at its bounds (T2f = 0.1 or 2 ms,
%      T2s = 5 or 20 ms, see T2estimate.m); exclBounds drops those from the fast/slow scatter
%    * the wrapped (aliased) medium signal adds a fast-decaying component to every 613 voxel,
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
cd(fullfile(fileparts(mfilename('fullpath')),'..','data','DWI_scan613'));   % release layout: run in the dataset folder
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','helper'));   % release layout: readT2Bounds.m
clc; clear; close all;

%% ---------------- 0. settings ----------------------------------------------------
DTI = pwd;                                          % .../DWI_scan613
T2D = fullfile(fileparts(pwd),'Na23_scan402-418');       % .../Na23_scan402-418
RAW = fullfile('..','bruker','20250801_201306_Mouse_23Naplus1_08012025_1_81');   % release layout
scanD = 613;   % proton DtiStandard (this folder)
scanT = 402;   % sodium scan that defines the Na23_scan402-418 grid

useRoll       = false;     % also roll about the tube axis so the organoid's off-axis direction matches
shiftMode     = 'organoid';% 'organoid': organoid centroids coincide (default)
                           % 'tube'    : tube axes coincide, the organoid only fixes the position along the axis
                           % 'both'    : half way between the two
extraShift_mm = [-0.1 0 0]; % manual tweak of the sodium data, mm along proton x(R-L) y(S-I) z(A-P); rerun 4-6
                           % [-0.1 0 0] = sodium maps 0.1 mm (1 proton voxel) to the LEFT in the coronal maps of 6b (-x = lower i).
                           % 2026-09-21: first set to -0.2 (central highest T2* long voxels on the high-ADC island of the core),
                           % then moved back by 0.1 mm because -0.2 looked too far in the overlay.  Step 4b prints the ADC under the
                           % central voxels for every shift from -0.4 to +0.2 mm (-0.2 and -0.1 give the same, highest, mean ADC)
T2set         = 'all';     % which sodium maps: 'all' = T2all / T2alls / T2allf (whole-matrix fit, mask ROIall, 2nd loop of T2estimate.m)
                           %                    'roi' = T2 / T2s / T2f (1st loop; since 2026-09-21 its mask is allave > 12000 = the
                           %                            medium, which leaves most of the organoid empty -> use 'all')
useSodiumROI  = false;     % false: the DWI ROI (roiD) alone defines the organoid area for ADC AND T2* (same area in both);
                           % true : additionally require the sodium ROI (Na23_scan402-418/ROI+orig), as before 2026-09-21
erodeROI      = 1;         % erode the proton ROI by n voxels (0.1 mm) before averaging ADC into sodium voxels
minFrac       = 0.5;       % sodium voxel must be covered >= this fraction by the (eroded) proton ROI
exclBounds    = true;      % fast / slow scatter: drop the sodium voxels whose T2f / T2s sits on its OWN fit bound
T2fBounds     = [];        % [] = read from T2estimate_scan402_418.m with readT2Bounds.m (the fit that belongs to T2set);
T2sBounds     = [];        %      or set by hand, e.g. [0.1 2] and [2 50] (ms)
mapMargin     = 3;         % proton-grid maps (6b): show the organoid ROI box (ROI+orig, all slices) + this many proton
                           % voxels (0.1 mm) around it; [] = full FOV

fprintf('proton folder : %s\nsodium folder : %s\n',DTI,T2D);

%% ---------------- 1. load the AFNI datasets and have a look ----------------------
[anatD,InfoD] = BrikLoad(fullfile(DTI,'allave+orig.BRIK'));      anatD = double(anatD);   % 80x64x64
roiD  = BrikLoad(fullfile(DTI,'ROI_shrink+orig.BRIK')) > 0;
adcD  = double(BrikLoad(fullfile(DTI,'ADC+orig.BRIK')));          % ADC in 1e-6 mm^2/s (myADCProt.m)
[anatT,InfoT] = BrikLoad(fullfile(T2D,'allave+orig.BRIK'));      anatT = double(anatT);
roiT  = BrikLoad(fullfile(T2D,'ROI+orig.BRIK')) > 0;
switch T2set
    case 'all', T2names = {'T2all','T2alls','T2allf'};  whichFit = 'last';    % whole matrix (ROIall)
    case 'roi', T2names = {'T2','T2s','T2f'};           whichFit = 'first';
    otherwise, error('T2set must be all or roi');
end
T2m   = double(BrikLoad(fullfile(T2D,[T2names{1} '+orig.BRIK'])));   % mono-exponential T2* (ms)
T2s   = double(BrikLoad(fullfile(T2D,[T2names{2} '+orig.BRIK'])));   % slow component (ms)
T2f   = double(BrikLoad(fullfile(T2D,[T2names{3} '+orig.BRIK'])));   % fast component (ms)
fprintf('sodium maps   : %s / %s / %s (%d voxels with a fit; %d of the %d sodium-ROI voxels)\n', T2names{:}, nnz(T2m>0), nnz(T2m>0 & roiT), nnz(roiT));
if isempty(T2fBounds) || isempty(T2sBounds)
    [T2fBounds, T2sBounds, bnd] = readT2Bounds(fullfile('..','..','scripts','T2estimate_scan402_418.m'), whichFit);
    fprintf('fit bounds read from T2estimate.m line %d (%s fit of %d): T2f [%g %g], T2s [%g %g] ms\n', bnd.line, whichFit, bnd.n, T2fBounds, T2sBounds);
end
ND = size(anatD); NT = size(anatT);
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
%  proton 613New : AFNI (i,j,k) = 2dseq (79-d3, d2, d1)  (proc613_to3d.sh: true labels + 3daxialize RSA)
%  sodium 402    : AFNI (i,j,k) = 2dseq (d1, d2, d3)      (generic to3d labels, no permutation)
PD = [0 0 1 0; 0 1 0 0; -1 0 0 ND(1)-1; 0 0 0 1];     % ijk -> d1 d2 d3 for the proton
PT = eye(4);                                            % ijk -> d1 d2 d3 for the sodium
[A_D,hdrD] = brukerIJK2subject(fullfile(RAW,num2str(scanD),'pdata','1','visu_pars'), PD);
[A_T,hdrT] = brukerIJK2subject(fullfile(RAW,num2str(scanT),'pdata','1','visu_pars'), PT);
describeFrame('proton 613', A_D, ND, hdrD);
describeFrame('sodium 402', A_T, NT, hdrT);
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
print(gcf,fullfile(DTI,'QC_tubefit_613New.jpg'),'-djpeg','-r150');

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
save(fullfile(DTI,'T2_to_613New_xform.mat'),'Q','A_D','A_T','M_T2ijk_to_Dijk','M_Dijk_to_T2ijk','dD','dT','cD','cT','tcD','tcT','useRoll','shiftMode','extraShift_mm');

%% 4b. optional check of extraShift_mm(1): ADC under the central highest T2* long voxels vs the left-right shift
%  The two sodium voxels with the highest T2* long at the sodium organoid centroid (centroid voxel + its neighbours) are
%  followed while the sodium data are moved along proton x in steps of one proton voxel (y, z of extraShift_mm kept).
%  ADC = box average of the DWI-ROI voxels that fall into the sodium voxel.  Only the sodium maps and the ROI centroid
%  choose the voxels, the ADC is not used for that.
Q0 = Q;  Q0(1:3,4) = t + A_D(1:3,1:3)*([0; extraShift_mm(2:3)']/resD);  M0 = A_D \ (Q0*A_T);      % alignment without the x tweak
[gi,gj,gk] = ndgrid(1:NT(1),1:NT(2),1:NT(3));
cand = find(sqrt((gi-cT_ijk(1)).^2 + (gj-cT_ijk(2)).^2 + (gk-cT_ijk(3)).^2) <= 1.2 & T2s > 0);
[~,o] = sort(T2s(cand),'descend');  cand = cand(o(1:min(2,numel(o))));  [c1,c2,c3] = ind2sub(NT,cand);
[I,J,K] = ndgrid(0:ND(1)-1,0:ND(2)-1,0:ND(3)-1);  P = [I(:)';J(:)';K(:)';ones(1,numel(I))];  clear I J K gi gj gk
fprintf('\nShift check (sodium organoid centroid = [%.1f %.1f %.1f]): mean ADC (um^2/s) of the DWI-ROI voxels under\n', cT_ijk);
for c = 1:numel(cand), fprintf('   voxel %d = sodium [%d %d %d], T2* long %.1f ms\n', c, c1(c), c2(c), c3(c), T2s(cand(c))); end
fprintf('   %8s', 'dx (mm)'); fprintf('   voxel %d', 1:numel(cand)); fprintf('     mean   proton i of voxel 1\n');
for dx = -0.4:0.1:0.2
    Sh = eye(4);  Sh(1,4) = dx/resD;  tr = round((Sh*M0) \ P) + 1;                  % proton voxel -> sodium voxel (1-based)
    a = nan(1,numel(cand));  ir = [NaN NaN];
    for c = 1:numel(cand)
        m = tr(1,:)' == c1(c) & tr(2,:)' == c2(c) & tr(3,:)' == c3(c) & roiD(:) & adcD(:) > 0;
        if any(m), a(c) = mean(adcD(m)); end
        if c == 1 && any(m), [pi1,~,~] = ind2sub(ND,find(m)); ir = [min(pi1) max(pi1)]; end
    end
    fprintf('   %+8.1f', dx); fprintf('   %7.0f', a); fprintf('  %7.0f   %d-%d%s\n', mean(a,'omitnan'), ir, ternaryStr(abs(dx-extraShift_mm(1))<1e-9,'   <- extraShift_mm',''));
end
clear P tr

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
writeOut(T2m_on ,InfoD,fullfile(DTI,'T2_on613New'));
writeOut(T2s_on ,InfoD,fullfile(DTI,'T2s_on613New'));
writeOut(T2f_on ,InfoD,fullfile(DTI,'T2f_on613New'));
writeOut(double(roiT_on),InfoD,fullfile(DTI,'ROI_T2_on613New'));
writeOut(anatT_on,InfoD,fullfile(DTI,'allave_T2_on613New'));

% 5b. proton ADC -> sodium grid (box average of the proton ROI voxels that fall in each sodium voxel)
roiDe = roiD; for e = 1:erodeROI, roiDe = imerode(roiDe, strel('sphere',1)); end
rv = double(roiDe(in)); rv = rv(:);  av = adcD(in); av = av(:);
cntAll = accumarray(lin(:),1,[prod(NT) 1]);                          % proton voxels per sodium voxel (~27)
cntRoi = accumarray(lin(:),rv,[prod(NT) 1]);                         % of which inside the eroded proton ROI
sumAdc = accumarray(lin(:),rv.*av,[prod(NT) 1]);
ADC_onT  = reshape(sumAdc./max(cntRoi,1),NT);  ADC_onT(cntRoi==0) = 0;
frac_onT = reshape(cntRoi./max(cntAll,1),NT);
writeOut(ADC_onT ,InfoT,fullfile(T2D,'ADC_613New_on402grid'));       % sodium grid -> next to the sodium allave
writeOut(frac_onT,InfoT,fullfile(T2D,'ROIfrac_613New_on402grid'));

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
print(gcf,fullfile(DTI,'QC_alignment_613New.jpg'),'-djpeg','-r150');

%% ---------------- 6. voxel-wise ADC vs T2* ------------------------------------------
% 6a. on the sodium grid (independent 0.3 mm voxels; each has its own T2* fit and the mean ADC of the
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
print(gcf,fullfile(DTI,'ADCvsT2_613New_sodiumgrid.jpg'),'-djpeg','-r150');
[si,sj,sk] = ind2sub(NT,find(sel));
writetable(table(si,sj,sk,ADC_onT(sel),T2f(sel),T2s(sel),T2m(sel),frac_onT(sel),atBound(sel), ...
           'VariableNames',{'i','j','k','ADC_1e-6mm2s','T2fast_ms','T2slow_ms','T2mono_ms','ROIcoverage','fitOnBound'}), ...
           fullfile(DTI,'ADCvsT2_613New_voxels.csv'));
fprintf('  wrote ADCvsT2_613New_voxels.csv\n');

%% 6b. on the proton grid (0.1 mm; T2* is nearest-neighbour, i.e. constant over 3x3x3 blocks) -> maps only
%  ONE mask for the three maps = the DWI ROI (roiD): ADC and T2* cover exactly the same area.  With the whole-matrix
%  T2* maps (T2set = 'all') every DWI-ROI voxel has a sodium fit, so nothing is lost; the count is printed.
selD = roiD & T2s_on>0 & T2f_on>0 & adcD>0;
kk=14;
fprintf('\nproton-grid maps, slice k=%d: DWI ROI %d voxels, with ADC and T2* %d (all slices: %d of %d)\n', kk, nnz(roiD(:,:,kk)), nnz(selD(:,:,kk)), nnz(selD), nnz(roiD));
roiBox = BrikLoad(fullfile(DTI,'ROI+orig.BRIK')) > 0;     % main organoid ROI: defines the field of view of the maps (mapMargin)
figure('Name','6 maps on proton grid','Position',[50 500 1500 400]);
subplot(1,3,1); showMap(anatD(:,:,kk), adcD(:,:,kk), selD(:,:,kk), 'ADC (um^2/s)',[100 300]); zoomToROI(roiBox,mapMargin); title(sprintf('ADC, coronal k=%d',kk));
subplot(1,3,2); showMap(anatD(:,:,kk), T2s_on(:,:,kk), selD(:,:,kk), 'T2* long (ms)',[0 10]); zoomToROI(roiBox,mapMargin); title('sodium T2* slow on proton grid');
subplot(1,3,3); showMap(anatD(:,:,kk), T2f_on(:,:,kk), selD(:,:,kk), 'T2* short (ms)',[0 1.5]); zoomToROI(roiBox,mapMargin); title('sodium T2* fast on proton grid');
title(sprintf('proton-grid pairs (n=%d, blocky T2*)',nnz(selD)),'FontSize',9);
set(findall(gcf,'Type','axes'),'Toolbar',[]);              % no axes-toolbar icons in a batch-mode export
print(gcf,fullfile(DTI,['ADCvsT2_613New_protongrid_slice',num2str(kk),'.jpg']),'-djpeg','-r150');

%% ---------------- 7. optional: did the sample move during the 613 scan? -------------
% organoid centre in each b-value volume (allave_t).  A drift of >1 voxel (0.1 mm) between frames
% means the per-voxel ADC fit mixes neighbouring voxels -> trust the eroded/box-averaged values more.
dat = double(BrikLoad(fullfile(DTI,'allave_t+orig.BRIK')));      % 80x64x64x7
dro = imdilate(roiD, strel('sphere',3));
fprintf('\nOrganoid centre per b-value volume (voxels, 1-based):\n');
for t = 1:size(dat,4)
    v = dat(:,:,:,t); thr = prctile(v(dro),60); m = v>thr & dro;
    [a,b,c] = ind2sub(ND,find(m)); fprintf('  frame %d : i=%.2f j=%.2f k=%.2f\n',t,mean(a),mean(b),mean(c));
end


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

function s = ternaryStr(c,a,b), if c, s = a; else, s = b; end, end

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
