% Fig. 4 add-on: R2 of the voxel-wise bi-exponential T2* fit within the ROI shown on the maps
% for organoid #3 (M1) and #4 (M3).
% Dots = individual ROI voxels. A box plot of the same voxels is written as R2bi_ROIbox_fig4.
% Needs ROI_T2bi_R2stats.mat in each folder, written by the last section of
% T2estimate_M1T2.m and T2estimate_M3T2.m (R2 = gof.rsquare from [result,gof]=fit(...)).
cd(fullfile(fileparts(mfilename('fullpath')),'..','..','data'));   % release layout: run in the dataset folder
clc;clear all;close all;
fs=16; % figure style (Final Figures.pptx): no grid, left + bottom axes only, thick black axes, black bold Arial
folders={'M1','M3'}; labels={'#3','#4'};
mu=zeros(1,2);sd=zeros(1,2);se=zeros(1,2);n=zeros(1,2);R2all=cell(1,2);
for f=1:2
    S=load([folders{f},'/ROI_T2bi_R2stats.mat']);
    R2all{f}=S.R2_ROI(:);
    n(f)=numel(S.R2_ROI); mu(f)=mean(S.R2_ROI); sd(f)=std(S.R2_ROI); se(f)=sd(f)/sqrt(n(f));
    fprintf('%s (%s, n=%d voxels), mean +/- std: T2*short=%.4f +/- %.4f ms, T2*long=%.4f +/- %.4f ms, R2=%.4f +/- %.4f (SEM=%.4f)\n',...
        labels{f},folders{f},n(f),mean(S.T2short_ROI),std(S.T2short_ROI),mean(S.T2long_ROI),std(S.T2long_ROI),mu(f),sd(f),se(f));
end
figure(Position=[100 100 340 340],Color='w');
bar(1:2,mu,0.6,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k');
hold on
for f=1:2 % one dot per ROI voxel
    swarmchart(f*ones(n(f),1),R2all{f},5,'k','filled','MarkerFaceAlpha',0.35,'MarkerEdgeColor','none','XJitterWidth',0.35);
end
errorbar(1:2,mu,se,'k','LineStyle','none','LineWidth',2,'CapSize',12);
xticks(1:2); xticklabels(labels); xlim([0.4 2.6]); ylim([0 1]); box off
ylabel('R^2 (bi-exponential fit)')
fontsize(gcf,fs,"points"); styleAx(gca)
print(gcf,[pwd,'/R2bi_ROIbar_fig4.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/R2bi_ROIbar_fig4.eps'],'-depsc','-r300');

% box plot of the same voxels (box = 25-75%, line = median, whiskers = 1.5 IQR, dots = voxels)
figure(Position=[100 100 340 340],Color='w');
xg=[ones(n(1),1);2*ones(n(2),1)]; yg=[R2all{1};R2all{2}];
boxchart(xg,yg,'BoxFaceColor',[0.75 0.75 0.75],'BoxFaceAlpha',1,'BoxEdgeColor','k','WhiskerLineColor','k','MarkerStyle','none','BoxWidth',0.6,'LineWidth',2);
hold on
swarmchart(xg,yg,10,'k','filled','MarkerFaceAlpha',0.35,'MarkerEdgeColor','none','XJitterWidth',0.35);
xticks(1:2); xticklabels(labels); xlim([0.4 2.6]); ylim([0.9 1]); box off
ylabel('Voxel-wise R^2')
fontsize(gcf,fs,"points"); styleAx(gca)
for f=1:2, fprintf('%s: R2 median = %.4f, 25-75%% = [%.4f %.4f]\n',labels{f},median(R2all{f}),prctile(R2all{f},[25 75])); end
print(gcf,[pwd,'/R2bi_ROIbox_fig4.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/R2bi_ROIbox_fig4.eps'],'-depsc','-r300');

%% local functions
function styleAx(ax)
    grid(ax,'off'); box(ax,'off');
    set(ax,'XColor','k','YColor','k','LineWidth',2,'FontName','Arial','FontWeight','bold');
    set([ax.XLabel ax.YLabel ax.Title],'Color','k','FontName','Arial','FontWeight','bold');
end
