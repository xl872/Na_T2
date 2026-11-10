cd(fullfile(fileparts(mfilename('fullpath')),'..','..','data','M3'));   % release layout: run in the dataset folder
clc;clear all;close all;
% figure style (Final Figures.pptx): no grid, left + bottom axes only, thick black axes, black bold Arial
fs=16; yLimPlot=[0 20]; yTickPlot=0:5:20;

time=[0.68 0.84 1 1.25 1.5 2 3 4 5 6.5 8 10 15 20 30 40 50 65 80];

% mask(mask<12000)=0;
% mask(mask>0)=1;
data=[];BG=[];
for run=[61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79] %[24 27 30 33 36 39 42 46 49 52 55 58]
    [V,info]=BrikLoad(['20250518_144207_Mouse_23Naplus1_05182025_02_1_66.',num2str(run),'.an200epi+orig.BRIK']);
    data0=squeeze(V);
    
    % loads a series of data with a number to string conversion
    
    %filter size of 3 (2*ceil(2*sigma)+1
    bg=data0([1:3,14:16],[1:3,14:16],:);
    BG=cat(2,BG,bg(:));
    data0=imgaussfilt3(data0,0.5);%-mean(bg,'all');
    %soothes out the voxels using 3D gaussian filter with 0.5 sigma for the
    data=cat(4,data,data0);
    
    % concatenates the data with the dimenson of time(4)
end



%% mask=1
mask=BrikLoad('ROI+orig.BRIK');
data_std=std(BG);
endn=length(time);
%amount of time points we want to use

data1=data(:,:,:,1:endn);
data1=data1.*repmat(mask,1,1,1,size(data1,4));
tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")./std(BG)';
% tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")';
% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0.1,0,0,5],...
   'Upper',[Inf,5, Inf,1,100],...
   'StartPoint',[tdata(1) 5 tdata(end) 0.5 5]);
ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
[result gof1]= fit(time(1:end)', tdata(1:end),ft);
result.b
result.e
%
result1=result;tdata1=tdata;

% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0],...
   'Upper',[Inf,50, Inf],...
   'StartPoint',[tdata(1) 5 tdata(end)]);
ft = fittype('c+a*exp(-x/b)','options',fo); 
[result,gof11] = fit(time(1:end)', tdata(1:end),ft);
result.b
result11=result;tdata11=tdata;
%

figure(Position=[100 100 476 300],Color='w');   % aspect of the slide frame (3.0 x 1.89 in)
hold on
te=linspace(time(1),time(end),500)';   % curves over the measured TEs, as in SNRfit_Fig4B.m (the short component rises steeply below the first TE)
plot(te,result11(te),'r','LineWidth',2)   % mono-exponential; fitting curves 2 pt, as on all supplementary figures
scatter(time(1:end)', tdata1(1:end),25,'k','filled')
plot(te,result1(te),'b','LineWidth',2)    % bi-exponential
xlim([0 time(end)]); xticks(0:20:80)

legend off                         % no legend (red = mono-, blue = bi-exponential; r2 printed below)
ylim(yLimPlot); yticks(yTickPlot)
ylabel('SNR')
xlabel('TE (ms)')
fontsize(gcf,fs,"points")
styleAx(gca)
fprintf('mono r2=%.4f  bi r2=%.4f | T2=%.3f ms, T2short=%.3f ms, T2long=%.3f ms\n',gof11.rsquare,gof1.rsquare,result11.b,result1.b,result1.e);
print(gcf,[pwd,'/aveT2_1&2fitted.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/aveT2_1&2fitted.eps'],'-depsc','-r300');
result1.d
%% mask=1
mask=BrikLoad('ROI_metrics+orig.BRIK');
data_std=std(BG);
endn=length(time);
%amount of time points we want to use

data1=data(:,:,:,1:endn);
data1=data1.*repmat(mask,1,1,1,size(data1,4));
tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")./std(BG)';
% tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")';
% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0.1,0,0,5],...
   'Upper',[Inf,5, Inf,1,100],...
   'StartPoint',[tdata(1) 0.1 tdata(end) 0.5 5]);
ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
[result gof1]= fit(time(1:end)', tdata(1:end),ft);
result.b
result.e
%
result1=result;tdata1=tdata;

% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0],...
   'Upper',[Inf,50, Inf],...
   'StartPoint',[tdata(1) 5 tdata(end)]);
ft = fittype('c+a*exp(-x/b)','options',fo); 
[result,gof11] = fit(time(1:end)', tdata(1:end),ft);
result.b
result11=result;tdata11=tdata;
%

figure(Position=[100 100 476 300],Color='w');   % aspect of the slide frame (3.0 x 1.89 in)
hold on
te=linspace(time(1),time(end),500)';   % curves over the measured TEs, as in SNRfit_Fig4B.m (the short component rises steeply below the first TE)
plot(te,result11(te),'r','LineWidth',2)   % mono-exponential; fitting curves 2 pt, as on all supplementary figures
scatter(time(1:end)', tdata1(1:end),25,'k','filled')
plot(te,result1(te),'b','LineWidth',2)    % bi-exponential
xlim([0 time(end)]); xticks(0:20:80)

legend off                         % no legend (red = mono-, blue = bi-exponential; r2 printed below)
ylim(yLimPlot); yticks(yTickPlot)
ylabel('SNR')
xlabel('TE (ms)')
fontsize(gcf,fs,"points")
styleAx(gca)
fprintf('mono r2=%.4f  bi r2=%.4f | T2=%.3f ms, T2short=%.3f ms, T2long=%.3f ms\n',gof11.rsquare,gof1.rsquare,result11.b,result1.b,result1.e);
print(gcf,[pwd,'/aveT2_1&2fitted_metrics.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/aveT2_1&2fitted_metrics.eps'],'-depsc','-r300');
result1.d
%% mean and std of the voxel-wise T2 maps within each ROI
% maps are T2/T2short/T2long+orig written by the voxel-wise T2estimate script
for roi={'ROI','ROI_metrics'}
    mask=BrikLoad([roi{1},'+orig.BRIK']);
    for map={'T2','T2short','T2long'}
        V=BrikLoad([map{1},'+orig.BRIK']);
        v=V(mask>0);
        fprintf('%s %s: mean=%.4f ms, std=%.4f ms, median=%.4f ms, 25-75%%=[%.4f %.4f] ms (n=%d voxels)\n',roi{1},map{1},mean(v),std(v),median(v),prctile(v,[25 75]),numel(v));
    end
end

%% local functions
function styleAx(ax)
    grid(ax,'off'); box(ax,'off');
    set(ax,'XColor','k','YColor','k','LineWidth',2,'FontName','Arial','FontWeight','bold');
    set([ax.XLabel ax.YLabel ax.Title],'Color','k','FontName','Arial','FontWeight','bold');
end
