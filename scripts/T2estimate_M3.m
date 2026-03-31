clc;clear all;close all;
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

%%
figure
imagesc(data(:,:,10,1)); axis xy
%%
figure
imagesc(mask(:,:,10)); axis xy
%%
figure;
subplot(2,1,1)
boxplot(BG);
xticklabels(num2str(time'))
ylabel('Image amplitude')
title('Backgroud value')
subplot(2,1,2)
scatter(time,std(BG),'filled');
ylabel('Backgroud STD')
xlabel('TE (s)')
xlim("tight")
fontsize(gcf,12,"points")
print(gcf,[pwd,'\Baseline.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\Baseline.eps'],'-depsc','-r300');
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
   'StartPoint',[tdata(1) 5 tdata(end) tdata(1) 5]);
ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
result = fit(time(1:end)', tdata(1:end),ft);
result.b
result.e
%
result1=result;tdata1=tdata;
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['T2short=',num2str(result.b),'ms ','T2long=',num2str(result.e),'ms'])
print(gcf,[pwd,'\AveT2_2estimateROI.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\AveT2_2estimateROI.eps'],'-depsc','-r300');
%
save('ROIT2estimate2.mat', 'time','tdata','result')

% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0],...
   'Upper',[Inf,50, Inf],...
   'StartPoint',[tdata(1) 5 tdata(end)]);
ft = fittype('c+a*exp(-x/b)','options',fo); 
result = fit(time(1:end)', tdata(1:end),ft);
result.b
result11=result;tdata11=tdata;
%
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['T2=',num2str(result.b),'ms'])
% legend('Location', 'southeast');
print(gcf,[pwd,'\AveT2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\AveT2estimate.eps'],'-depsc','-r300');


%%
mask=BrikLoad('formask+orig.BRIK');
mask(mask<10000)=0;
mask(mask>0)=1;
maskCor=mask;
% maskCor(:,7:end,:)=0;
maskCSF=BrikLoad('ROI+orig.BRIK');
data=data(:,:,:,:);
data_std=std(BG);
T2=[];T2f=[];T2s=[];Tdata1=[];Tdata2=[];
for i=1:size(data,1)
    for j=1:size(data,2)
       for k=1:size(data,3)
            if mask(i,j,k)~=0
                tdata=squeeze(data(i,j,k,:))./data_std';
               if maskCor(i,j,k)==1
                    Tdata1=[Tdata1 tdata];
                elseif maskCSF(i,j,k)==1
                    Tdata2=[Tdata2 tdata];
                end
                fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0,0,5],...
                   'Upper',[Inf,5, Inf,1,100],...
                   'StartPoint',[tdata(1) 1 tdata(end) tdata(1) 5]);
                ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
                result = fit(time(1:end)', tdata(1:end),ft);
                T2f(i,j,k)=min(result.b,result.e);
                T2s(i,j,k)=max(result.b,result.e);

              
                fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[Inf,100, Inf],...
                   'StartPoint',[tdata(1) 5 tdata(end)]);
                ft = fittype('c+a*exp(-x/b)','options',fo); 
                result = fit(time(1:end)', tdata(1:end),ft);
                T2(i,j,k)=result.b;
            else
                T2f(i,j,k)=0;
                T2s(i,j,k)=0;
                T2(i,j,k)=0;
            end

       end
    end
end

%%
info.BRICK_TYPES=3;
opt.OverWrite='y';
opt.prefix='T2';
WriteBrik(T2,info,opt);
opt.prefix='T2short';
WriteBrik(T2f,info,opt);
opt.prefix='T2long';
WriteBrik(T2s,info,opt);

opt.prefix='T2masked';
WriteBrik(T2.*maskCSF,info,opt);
opt.prefix='T2shortmasked';
WriteBrik(T2f.*maskCSF,info,opt);
opt.prefix='T2longmasked';
WriteBrik(T2s.*maskCSF,info,opt);
%%
%%
% anat=squeeze(BrikLoad('20250518_144207_Mouse_23Naplus1_05182025_02_1_66.80.an200epi.blur.resample+orig.BRIK'));
anat=squeeze(BrikLoad('formask+orig.BRIK'));

% maskt=BrikLoad('COPY_20250101Na.13.an200epi.resamp+orig');
% maskt(:,80:end,:)=0;
% mask=BrikLoad('ROI+orig.BRIK');

mask=BrikLoad('ROI+orig.BRIK');
% mask(mask<12000)=0;
% mask(mask>0)=1;

maskt=mask;

for sl=9:11%:size(data,3)
figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
B=mean(anat(:,:,sl),3)';
B=repmat(B,1,1,3)/max(B(:));
hB = image(B);%axis xy
axis image off;
hold on;
F=mean(T2f(:,:,sl),3)';
% F=imresize(F,[size(B,1),size(B,2)]);
hF = imagesc(F);%axis xy
colormap('winter')
clim([0 0.5]);
%
alphadata =1*min(maskt(:,:,sl),[],3)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2*short (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2shortestimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2shortestimate',num2str(sl),'.eps'],'-depsc','-r300');


figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
% B=anat(:,:,sl)';
% B=repmat(B,1,1,3)/max(B(:))*1.5;
hB = image(B);%axis xy
axis image off;
hold on;
F=mean(T2s(:,:,sl),3)';
% F=imresize(F,[size(B,1),size(B,2)],"bilinear");
hF = imagesc(F);%axis xy
colormap('autumn')
clim([5 20]);
% alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2*long (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2longestimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2longestimate',num2str(sl),'.eps'],'-depsc','-r300');

figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
% B=anat(:,:,sl)';
% B=repmat(B,1,1,3)/max(B(:))*1.5;
hB = image(B);%axis xy
axis image off;
hold on;
F=mean(T2(:,:,sl),3)';
% F=imresize(F,[size(B,1),size(B,2)],"bilinear");
hF = imagesc(F);%axis xy
colormap('autumn')
clim([5 20]);
% alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2* (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2estimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2estimate',num2str(sl),'.eps'],'-depsc','-r300');

figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
% B=anat(:,:,sl)';
% B=repmat(B,1,1,3)/max(B(:))*1.5;
hB = image(B);%axis xy
axis image off;
hold on;
F=mean(T2(:,:,sl),3)';
% F=imresize(F,[size(B,1),size(B,2)],"bilinear");
hF = imagesc(F);%axis xy
colormap('autumn')
clim([5 20]);
% alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',0);
h = colorbar;
ylabel(h, 'T2* (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2estimate',num2str(sl),'BG.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2estimate',num2str(sl),'BG.eps'],'-depsc','-r300');


end




%% mask=1
mask=BrikLoad('ROIl+orig.BRIK');
% mask(:,:,[11])=0;
data_std=std(BG);
endn=length(time);
%amount of time points we want to use

data1=data(:,:,:,1:endn);
data1=data1.*repmat(mask,1,1,1,size(data1,4));
tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")./std(BG)';

% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0.1,0,0,5],...
                   'Upper',[Inf,5, Inf,1,100],...
                   'StartPoint',[tdata(1) 5 tdata(end) tdata(1) 5]);
                ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
                result = fit(time(1:end)', tdata(1:end),ft);
                
result.b
result.e
%
result1=result;tdata1=tdata;
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['T2short=',num2str(result.b),'ms ','T2long=',num2str(result.e),'ms'])
print(gcf,[pwd,'\AveT2_2estimateROIl.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\AveT2_2estimateROIl.eps'],'-depsc','-r300');
%
mask=BrikLoad('ROIh+orig.BRIK');

data_std=std(BG);
endn=length(time);
%amount of time points we want to use

data1=data(:,:,:,1:endn);
data1=data1.*repmat(mask,1,1,1,size(data1,4));
tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")./std(BG)';

% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0.1,0,0,5],...
                   'Upper',[Inf,5, Inf,1,100],...
                   'StartPoint',[tdata(1) 5 tdata(end) tdata(1) 5]);
                ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
result = fit(time(1:end)', tdata(1:end),ft);
result.b
result.e
%
result1=result;tdata1=tdata;
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['T2short=',num2str(result.b),'ms ','T2long=',num2str(result.e),'ms'])
print(gcf,[pwd,'\AveT2_2estimateROIh.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\AveT2_2estimateROIh.eps'],'-depsc','-r300');