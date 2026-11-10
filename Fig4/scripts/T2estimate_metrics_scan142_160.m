% T2estimate_metrics_scan142_160.m - 20250801 session, runs 142-160
% Same process as page 2 of Organoid_supp_revised pptx (T2estimate_M3T2.m + T2estimate_metrics.m + R2mono_table.m):
%  1) ROI-mean SNR vs TE, mono- (red) and bi-exponential (blue) fit, for ROI+orig ("Organoid") and ROI_metrics+orig ("Metrics")
%       -> aveT2_1&2fitted.jpg/.eps (Organoid), aveT2_1&2fitted_metrics.jpg/.eps (Metrics), aveT2fit_ROIs.csv
%  2) voxel-wise mono- and bi-exponential T2* fits (fit settings of T2estimate_M3T2.m and R2mono_table.m)
%       -> map_T2mono, map_T2short, map_T2long, map_R2mono, map_R2bi +orig
%          (own prefixes: T2+orig/T2f/T2s/T2R2bi/T2R2mono of T2estimate.m are used by DWI_scan175 and are not touched)
%  3) mean, std and CV (= std/mean) of the voxel-wise maps within each ROI -> T2_ROI_stats.csv
%  4) T2* (mono) map of 3 slices with the outlines of the real ROI / ROI_metrics voxels -> T2map_ROIcontour.png/.jpg/.eps
cd(fullfile(fileparts(mfilename('fullpath')),'..','data','Na23_scan142-160'));   % release layout: run in the dataset folder
clc;clear all;close all;
addpath(genpath(fullfile(getenv('HOME'),'MATLAB Add-Ons','Toolboxes','afni')));   % BrikLoad / WriteBrik
fs=16;                    % fit figures: style of T2estimate_metrics.m (Final Figures.pptx)
pre='20250801_201306_Mouse_23Naplus1_08012025_1_81.';
runs=142:160;
time=[0.68 0.84 1 1.25 1.5 2 3 4 5 6.5 8 10 15 20 30 40 50 65 80];   % PVM_EchoTime (ms) of runs 142-160
sls=7:9;                  % slices shown on the map (ROI: slices 7-10, ROI_metrics: slices 6-8)
vox=0.4;                  % in-plane voxel size (mm), for the 2 mm scale bar
thr=10000;                % allave threshold of the whole-organoid mask shown on the map (as formask>=10000 for M1/M3)
clT2=[5 20];              % colour range of the T2* map (ms), as page 2
refit=true;               % false: reload the voxel-wise maps written by an earlier run
roiFile={'ROI','ROI_metrics'}; roiName={'Organoid','Metrics'};
roiCol=[25 107 36; 160 43 147]/255;   % pptx accent3 (green, "Organoid") / accent5 (purple, "Metrics")

data=[];BG=[];
for run=runs
    data0=squeeze(BrikLoad([pre,num2str(run),'.an200epi+orig.BRIK']));
    bg=data0([1:3,14:16],[1:3,14:16],:);    % background corners
    BG=cat(2,BG,bg(:));
    data0=imgaussfilt3(data0,0.5);           % 3D gaussian filter, sigma 0.5 voxel
    data=cat(4,data,data0);
end
data_std=std(BG);                            % noise per TE: SNR = signal / data_std

%% 1) ROI-mean SNR decay with mono- and bi-exponential fit (as T2estimate_metrics.m)
for r=1:2
    mask=BrikLoad([roiFile{r},'+orig.BRIK'])>0;
    tdata=squeeze(sum(data.*mask,[1,2,3]))/nnz(mask)./data_std';
    % bi-exponential: c + a*((1-d)*exp(-TE/T2short) + d*exp(-TE/T2long)); T2short start 5 ms (Organoid) / 0.1 ms (Metrics)
    b0=[5 0.1];
    fo = fitoptions('Method','NonlinearLeastSquares',...
       'Lower',[0,0.1,0,0,5],...
       'Upper',[Inf,5, Inf,1,100],...
       'StartPoint',[tdata(1) b0(r) tdata(end) 0.5 5]);
    ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo);
    [fbi,gbi,obi] = fit(time', tdata,ft);   % obi.exitflag 0 = stopped at the default 600 function evaluations (as the page-2 M3 fits)
    % mono-exponential: c + a*exp(-TE/T2)
    fo = fitoptions('Method','NonlinearLeastSquares',...
       'Lower',[0,0,0],...
       'Upper',[Inf,50, Inf],...
       'StartPoint',[tdata(1) 5 tdata(end)]);
    ft = fittype('c+a*exp(-x/b)','options',fo);
    [fmono,gmono] = fit(time', tdata,ft);
    ave(r)=struct('tdata',tdata,'fbi',fbi,'gbi',gbi,'fmono',fmono,'gmono',gmono,'exitflag',obi.exitflag); %#ok<SAGROW>
    fprintf('%-8s (n=%d) mono r2=%.4f  bi r2=%.4f | T2=%.3f ms, T2short=%.3f ms, T2long=%.3f ms, long fraction d=%.3f (bi exitflag %d)\n',...
        roiName{r},nnz(mask),gmono.rsquare,gbi.rsquare,fmono.b,fbi.b,fbi.e,fbi.d,obi.exitflag);
end
yTop=max(20,5*ceil(max([ave.tdata],[],'all')/5));   % same SNR axis for both ROIs (page 2: 0-20)
outFile={'aveT2_1&2fitted','aveT2_1&2fitted_metrics'};
for r=1:2
    figure(Position=[100 100 476 300],Color='w');   % aspect of the slide frame (3.0 x 1.89 in)
    hold on
    te=linspace(min(time),max(time),500)';
    plot(te,ave(r).fmono(te),'r','LineWidth',2)      % mono-exponential
    scatter(time', ave(r).tdata,25,'k','filled')
    plot(te,ave(r).fbi(te),'b','LineWidth',2)        % bi-exponential
    xlim([0 max(time)]); xticks(0:20:80)
    ylim([0 yTop]); yticks(0:5:yTop)
    ylabel('SNR')
    xlabel('TE (ms)')
    fontsize(gcf,fs,"points")
    styleAx(gca)
    print(gcf,fullfile(pwd,[outFile{r},'.jpg']),'-djpeg','-r300');
    print(gcf,fullfile(pwd,[outFile{r},'.eps']),'-depsc','-r300');
end
Ta=table(roiName',arrayfun(@(a) a.fmono.b,ave)',arrayfun(@(a) a.gmono.rsquare,ave)',arrayfun(@(a) a.fbi.b,ave)',...
    arrayfun(@(a) a.fbi.e,ave)',arrayfun(@(a) a.fbi.d,ave)',arrayfun(@(a) a.gbi.rsquare,ave)',[ave.exitflag]',...
    'VariableNames',{'ROI','T2mono_ms','R2mono','T2short_ms','T2long_ms','long_fraction','R2bi','bi_exitflag'});
% T2short and long_fraction are not determined when the short component has decayed before the first TE (0.68 ms):
% check bi_exitflag and whether T2short / long_fraction sit at a bound before quoting them.
writetable(Ta,'aveT2fit_ROIs.csv');

%% 2) voxel-wise T2* fits (as T2estimate_M3T2.m; R2mono as R2mono_table.m)
[anat,info]=BrikLoad('allave+orig.BRIK');
roi={BrikLoad('ROI+orig.BRIK')>0, BrikLoad('ROI_metrics+orig.BRIK')>0};
disk=anat>=thr;                          % whole organoid (coloured part of the map)
fitmask=disk|roi{1}|roi{2};              % every ROI voxel is fitted, also if it is below thr
mapName={'map_T2mono','map_T2short','map_T2long','map_R2mono','map_R2bi'};
if refit || ~all(cellfun(@(f) exist([f,'+orig.HEAD'],'file')==2,mapName))
    [T2,T2f,T2s,R2mono,R2bi]=voxelFit(data,data_std,time,fitmask);
    info.BRICK_TYPES=3; opt.OverWrite='y';
    M={T2,T2f,T2s,R2mono,R2bi};
    for m=1:5
        opt.prefix=mapName{m};
        WriteBrik(M{m},info,opt);
    end
else
    T2=BrikLoad('map_T2mono+orig.BRIK'); T2f=BrikLoad('map_T2short+orig.BRIK'); T2s=BrikLoad('map_T2long+orig.BRIK');
    R2mono=BrikLoad('map_R2mono+orig.BRIK'); R2bi=BrikLoad('map_R2bi+orig.BRIK');
end

%% 3) mean, std and CV (= std/mean) within each ROI, rows as the page-2 table
qName={'T2 mono','T2-mono R2','T2 long','T2 short','T2-bi R2'};
qMap={T2,R2mono,T2s,T2f,R2bi};
T=roiStats(qMap,qName,roi,roiName);
writetable(T,'T2_ROI_stats.csv');
printTable(T,roiName,'Supp. Tables 4-5 #5');

%% 4) T2* map with the outlines of the real ROI / ROI_metrics voxels
sl2d=@(V,s) V(:,:,s)';                   % display orientation of T2estimate.m (rows = dim 2, columns = dim 1)
U={};F={};A={};R={};
for s=1:numel(sls)
    U{s}=sl2d(anat,sls(s)); F{s}=sl2d(T2,sls(s)); A{s}=sl2d(disk,sls(s));
    R{s,1}=sl2d(roi{1},sls(s)); R{s,2}=sl2d(roi{2},sls(s));
end
mapMontage(U,F,A,R,roiCol,clT2,vox,'T2* (ms)','T2map_ROIcontour');


%% local functions
function [T2,T2f,T2s,R2mono,R2bi]=voxelFit(data,data_std,time,fitmask)
% voxel-wise bi-exponential (settings of T2estimate_M3T2.m) and mono-exponential fit (settings of R2mono_table.m)
sz=size(fitmask);
T2=zeros(sz); T2f=zeros(sz); T2s=zeros(sz); R2mono=zeros(sz); R2bi=zeros(sz);
D=reshape(data,[],numel(time));
idxList=find(fitmask)';
fprintf('voxel-wise fits of %d voxels ...\n',numel(idxList));
for idx=idxList
    tdata=D(idx,:)'./data_std';
    fo = fitoptions('Method','NonlinearLeastSquares',...
       'Lower',[0,0,0,0,5],...
       'Upper',[Inf,5, Inf,1,100],...
       'StartPoint',[tdata(1) 1 tdata(end) tdata(1) 5]);
    ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo);
    [result,gof] = fit(time', tdata,ft);
    T2f(idx)=min(result.b,result.e);
    T2s(idx)=max(result.b,result.e);
    R2bi(idx)=gof.rsquare;
    fo = fitoptions('Method','NonlinearLeastSquares',...
       'Lower',[0,0,0],...
       'Upper',[Inf,100, Inf],...
       'StartPoint',[tdata(1) 5 tdata(end)]);
    ft = fittype('c+a*exp(-x/b)','options',fo);
    [result,gof] = fit(time', tdata,ft);
    T2(idx)=result.b;
    R2mono(idx)=gof.rsquare;
end
end

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

function styleAx(ax)
    grid(ax,'off'); box(ax,'off');
    set(ax,'XColor','k','YColor','k','LineWidth',2,'FontName','Arial','FontWeight','bold');
    set([ax.XLabel ax.YLabel ax.Title],'Color','k','FontName','Arial','FontWeight','bold');
end
