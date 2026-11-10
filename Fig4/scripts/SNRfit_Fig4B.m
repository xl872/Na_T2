% Fig. 4B remake (first made in Python): mean SNR of the organoid ROI vs TE with the
% bi-exponential fit c + a*((1-d)*exp(-TE/T2short) + d*exp(-TE/T2long)), for organoid #3 (M1) and #4 (M3).
% Uses ROIT2estimate2.mat (time, tdata = ROI-mean SNR, result = the fit) written by T2estimate_M1T2.m and
% T2estimate_M3T2.m - the same data and fit as the "Organoid" plots of T2estimate_metrics.m.
% Style of the supplementary figures: no grid, left + bottom axes only, thick black axes, black bold Arial, 2-pt curves.
% Writes SNRfit_Fig4B_3.jpg/.eps and SNRfit_Fig4B_4.jpg/.eps in the data folder (the script changes to it itself).
cd(fullfile(fileparts(mfilename('fullpath')),'..','..','data'));   % release layout: run in the dataset folder
clc;clear;close all;
folders={'M1','M3'}; labels={'#3','#4'};
fs=16; yLimPlot=[0 20]; yTickPlot=0:5:20;
dotCol=[0 0 0]; fitCol=[0 0 1];   % black dots, blue bi-exponential curve (as the page-2 plots)
for f=1:2
    S=load([folders{f},'/ROIT2estimate2.mat']);
    te=linspace(S.time(1),S.time(end),500)';   % curve over the measured TEs (the short component rises steeply below the first TE)
    figure(Position=[100 100 476 330],Color='w'); hold on   % page-2 plot size (476 x 300) + room for the title
    plot(te,S.result(te),'-','Color',fitCol,'LineWidth',2);
    scatter(S.time,S.tdata,25,dotCol,'filled');
    xlim([0 S.time(end)]); ylim(yLimPlot); yticks(yTickPlot); xticks(0:20:80);
    xlabel('TE (ms)'); ylabel('SNR'); title(['Organoid ',labels{f}]);
    fontsize(gcf,fs,"points"); styleAx(gca);
    fprintf('Organoid %s (%s): T2*short = %.3f ms, T2*long = %.3f ms, long fraction d = %.3f\n',labels{f},folders{f},S.result.b,S.result.e,S.result.d);
    print(gcf,[pwd,'/SNRfit_Fig4B_',labels{f}(2:end),'.jpg'],'-djpeg','-r300');
    print(gcf,[pwd,'/SNRfit_Fig4B_',labels{f}(2:end),'.eps'],'-depsc','-r300');
end

%% local functions
function styleAx(ax)
    grid(ax,'off'); box(ax,'off');
    set(ax,'XColor','k','YColor','k','LineWidth',2,'FontName','Arial','FontWeight','bold');
    set([ax.XLabel ax.YLabel ax.Title],'Color','k','FontName','Arial','FontWeight','bold');
end
