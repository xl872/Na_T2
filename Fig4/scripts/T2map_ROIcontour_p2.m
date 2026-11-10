% T2map_ROIcontour_p2.m - page 2 of Organoid_supp_revised pptx, organoid #3 (M1) and #4 (M3):
%  1) mean, std and CV (= std/mean) of the voxel-wise maps within ROI+orig ("Organoid") and ROI_metrics+orig ("Metrics")
%       maps written by T2estimate_M1T2.m, T2estimate_M3T2.m (T2, T2short, T2long, R2bi) and R2mono_table.m (R2mono)
%       -> <dir>/T2_ROI_stats.csv and T2_ROI_stats_p2.csv (both organoids)
%  2) T2* (mono) map of the 3 slices of page 2 with the outlines of the real ROI / ROI_metrics voxels
%       (replaces the circles drawn in the pptx) -> <dir>/T2map_ROIcontour.png/.jpg/.eps
cd(fullfile(fileparts(mfilename('fullpath')),'..','..','data'));   % release layout: run in the dataset folder
clc;clear all;close all;
addpath(genpath(fullfile(getenv('HOME'),'MATLAB Add-Ons','Toolboxes','afni')));   % BrikLoad
base=pwd;
S(1).dir='M1'; S(1).name='#3'; S(1).sls=12:14; S(1).sdim=1; S(1).vox=0.3;   % slices along dim 1 (T2estimate_M1T2.m)
S(2).dir='M3'; S(2).name='#4'; S(2).sls=9:11;  S(2).sdim=3; S(2).vox=0.4;   % slices along dim 3 (T2estimate_M3T2.m)
thr=10000;                % formask threshold of the voxel-wise fits = coloured part of the map (T2allestimate*.jpg)
clT2=[5 20];              % colour range of the T2* map (ms), as page 2
roiFile={'ROI','ROI_metrics'}; roiName={'Organoid','Metrics'};
roiCol=[25 107 36; 160 43 147]/255;   % pptx accent3 (green, "Organoid") / accent5 (purple, "Metrics")
qName={'T2 mono','T2-mono R2','T2 long','T2 short','T2-bi R2'};           % rows of the page-2 table
qFile={'T2','R2mono','T2long','T2short','R2bi'};
Tall=table();
for s=1:2
    cd(fullfile(base,S(s).dir));
    anat=BrikLoad('formask+orig.BRIK');
    disk=anat>=thr;
    roi={BrikLoad('ROI+orig.BRIK')>0, BrikLoad('ROI_metrics+orig.BRIK')>0};
    qMap=cellfun(@(f) BrikLoad([f,'+orig.BRIK']),qFile,'UniformOutput',false);
    % 1) statistics
    T=roiStats(qMap,qName,roi,roiName);
    writetable(T,'T2_ROI_stats.csv');
    printTable(T,roiName,S(s).name);
    Tall=[Tall; [table(repmat({S(s).name},height(T),1),repmat({S(s).dir},height(T),1),'VariableNames',{'organoid','dir'}) T]]; %#ok<AGROW>
    % 2) T2* map with the ROI outlines
    if S(s).sdim==1, sl2d=@(V,k) squeeze(V(k,:,:))'; else, sl2d=@(V,k) V(:,:,k)'; end   % display orientation of the T2estimate scripts
    U={};F={};A={};R={};
    for k=1:numel(S(s).sls)
        sl=S(s).sls(k);
        U{k}=sl2d(anat,sl); F{k}=sl2d(qMap{1},sl); A{k}=sl2d(disk,sl);
        R{k,1}=sl2d(roi{1},sl); R{k,2}=sl2d(roi{2},sl);
    end
    mapMontage(U,F,A,R,roiCol,clT2,S(s).vox,'T2* (ms)','T2map_ROIcontour');
end
cd(base);
writetable(Tall,'T2_ROI_stats_p2.csv');


%% local functions
function T=roiStats(qMap,qName,roi,roiName)
% n, mean, std (n-1), CV = std/mean, median and quartiles of each map within each ROI
rows={};
for r=1:numel(roi)
    for q=1:numel(qMap)
        x=double(qMap{q}(roi{r}));
        rows(end+1,:)={roiName{r},qName{q},numel(x),mean(x),std(x),std(x)/mean(x),median(x),prctile(x,25),prctile(x,75)}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'ROI','quantity','n_voxels','mean','std','CV','median','p25','p75'});
end

function printTable(T,roiName,label)
% same layout as the table on the slide: quantity | Organoid mean std CV | Metrics mean std CV
fprintf('\n%-11s| %-8s %24s | %-8s %24s\n',label,roiName{1},'',roiName{2},'');
fprintf('%-11s| %8s %8s %8s | %8s %8s %8s\n','T2* (ms)','Mean','Std','CV','Mean','Std','CV');
q=unique(T.quantity,'stable');
for i=1:numel(q)
    a=T(strcmp(T.ROI,roiName{1})&strcmp(T.quantity,q{i}),:); b=T(strcmp(T.ROI,roiName{2})&strcmp(T.quantity,q{i}),:);
    fprintf('%-11s| %8.3f %8.3f %7.1f%% | %8.3f %8.3f %7.1f%%\n',q{i},a.mean,a.std,100*a.CV,b.mean,b.std,100*b.CV);
end
fprintf('(n = %d / %d voxels)\n',a.n_voxels,b.n_voxels);
end

function mapMontage(U,F,A,R,col,cl,vox,cbLabel,fname)
% slices side by side: grey underlay U, T2* map F (autumn) where A, outline of the voxels of ROI R{s,r} in colour col(r,:),
% 2 mm scale bar on the first slice, colour bar on the right. Sized for the slide (about 5.6 x 1.85 in).
ns=numel(U); [ny,nx]=size(U{1});
pw=1.62; ph=pw*ny/nx; gap=0.03; x0=0.02; y0=0.1;           % panel size, spacing and margins (inches)
cbx=x0+ns*pw+(ns-1)*gap+0.08; cbw=0.11;
fig=figure(Units='inches',Position=[1 1 cbx+cbw+0.62 ph+2*y0],Color='w',InvertHardcopy='off');
for s=1:ns
    ax=axes(fig,Units='inches',Position=[x0+(s-1)*(pw+gap) y0 pw ph]);
    B=double(U{s}); B=B/max(B(:));
    image(ax,repmat(B,1,1,3)); hold(ax,'on'); axis(ax,'image','off');
    hF=imagesc(ax,F{s}); set(hF,'AlphaData',double(A{s}));
    colormap(ax,autumn); clim(ax,cl);
    for r=1:size(R,2)
        roiOutline(ax,R{s,r},col(r,:),2.5);
    end
    if s==1   % 2 mm scale bar, bottom left
        xb=1.2; yb=ny-0.6;
        plot(ax,[xb xb+2/vox]-0.5,[yb yb],'w','LineWidth',3);
        text(ax,xb-0.5,yb-0.9,'2mm','Color','w','FontName','Arial','FontSize',10,'FontWeight','bold','VerticalAlignment','bottom');
    end
end
axc=axes(fig,Units='inches',Position=[cbx y0 cbw ph]);    % colour bar: gradient with ticks and label on the right
imagesc(axc,[0 1],cl,linspace(cl(1),cl(2),256)'); colormap(axc,autumn); clim(axc,cl);
set(axc,'YDir','normal','XTick',[],'XColor','none','YAxisLocation','right','YTick',cl(1):5:cl(2),'TickLength',[0 0],...
    'Box','off','YColor','k','LineWidth',0.5,'FontName','Arial','FontSize',11);
ylabel(axc,cbLabel,'FontName','Arial','FontSize',11,'Interpreter','none');
print(fig,fullfile(pwd,[fname,'.png']),'-dpng','-r300');
print(fig,fullfile(pwd,[fname,'.jpg']),'-djpeg','-r300');
print(fig,fullfile(pwd,[fname,'.eps']),'-depsc','-r300');
end

function roiOutline(ax,M,col,lw)
% outline of the voxels of the 2-D mask M (voxel (i,j) of the image spans j+-0.5, i+-0.5), drawn just inside the
% voxel edges (inset 0.22 voxel ~ half the line width) so that the outlines of two touching ROIs lie side by side
if ~any(M(:)), return; end
k=50; r=11;                             % trace the mask upsampled k times and eroded by r/k = 0.22 voxel
P=bwboundaries(imerode(kron(double(M),ones(k)),ones(2*r+1)));
for b=1:numel(P)
    plot(ax,(P{b}(:,2)-0.5)/k+0.5,(P{b}(:,1)-0.5)/k+0.5,'Color',col,'LineWidth',lw);
end
end
