cd(fullfile(fileparts(mfilename('fullpath')),'..','data','Na23_scan142-160'));   % release layout: run in the dataset folder
clc;clear all;close all;
time=[0.68 0.84 1 1.25 1.5 2 3 4 5 6.5 8 10 15 20 30 40 50 65 80];
mask=BrikLoad('ROI+orig.BRIK');
% mask(mask<12000)=0;
%% 
%% 
% mask(mask>0)=1;
data=[];BG=[];
for run=[142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160] %[24 27 30 33 36 39 42 46 49 52 55 58]
    data0=squeeze(BrikLoad(['20250801_201306_Mouse_23Naplus1_08012025_1_81.',num2str(run),'.an200epi+orig.BRIK']));
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
imagesc(data(:,:,8,7)); axis xy
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
print(gcf,[pwd,'/Baseline.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/Baseline.eps'],'-depsc','-r300');
%% mask=1
% mask=BrikLoad('myMASK+orig.BRIK');
% csf=BrikLoad('CSF+orig.BRIK');
% mask(:,7:end,:)=0;
% mask(csf==1)=0;
data_std=std(BG);
endn=length(time);
%amount of time points we want to use

data1=data(:,:,:,1:endn);
data1=data1.*repmat(mask,1,1,1,size(data1,4));
tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")./std(BG)';
% tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")';
% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0,0,5],...
   'Upper',[Inf,1.5, Inf,1,30],...
   'StartPoint',[tdata(1) 0.1 tdata(end) tdata(1) 5]);
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
print(gcf,[pwd,'/AveT2_2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/AveT2_2estimate.eps'],'-depsc','-r300');
%

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
print(gcf,[pwd,'/AveT2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/AveT2estimate.eps'],'-depsc','-r300');


%%
% 
% figure(Position=[100 100 400 300]);
% plot(result11,'r')
% hold on
% scatter(time(1:end)', tdata1(1:end),'r','.')
% % 
% % plot(result22,'b')
% % scatter(time(1:end)', tdata2(1:end),'b','.')
% 
% legend({'Cortex fitted','Cortex data','CSF fitted','CSF data'},'Box','off')
% fontsize(gcf,10,"points")
% ylabel('SNR')
% xlabel('TE (ms)')
% title({['CortexT2=',num2str(result11.b),'ms '] ...
%     ['CSFT2=',num2str(result22.b),'ms']})
% print(gcf,[pwd,'/aveT2_CSF&COTestimate1.jpg'],'-djpeg','-r300');
% print(gcf,[pwd,'/aveT2_CSF&COTestimate1.eps'],'-depsc','-r300');
%%
mask=BrikLoad('allave+orig.BRIK'); 
mask(mask<12000)=0;
mask(mask>0)=1;
maskCor=mask;
% maskCor(:,7:end,:)=0;
maskCSF=BrikLoad('ROI+orig.BRIK');
data=data(:,:,:,:);
T2=[];T2f=[];T2s=[];R2bi=[];R2mono=[];Tdata1=[];Tdata2=[];   % R2bi / R2mono = R^2 of the two fits per voxel
for i=1:size(data,1)
    for j=1:size(data,2)
       for k=1:size(data,3)
            if mask(i,j,k)~=0
                % fo = fitoptions('Method','NonlinearLeastSquares',...
                %    'Lower',[0,0,0,0,0],...
                %    'Upper',[Inf,100, Inf,Inf,100],...
                %    'StartPoint',[tdata(1)/2 0.01 mdata/2 tdata(1)/100 10]);
                % ft = fittype('c+a*exp(-x/b)+d*exp(-x/e)','options',fo); 
                % % contains addition to the function to calculate T2
                % result = fit(time', squeeze(data(i,j,k,:)),ft);
                % T2f(i,j,k)=min(result.b,result.e);
                % T2s(i,j,k)=max(result.b,result.e);
                tdata=squeeze(data(i,j,k,:))./data_std';
               if maskCor(i,j,k)==1
                    Tdata1=[Tdata1 tdata];
                elseif maskCSF(i,j,k)==1
                    Tdata2=[Tdata2 tdata];
                end
                fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0,0,5],...
                   'Upper',[Inf,5, Inf,1,50],...
                   'StartPoint',[tdata(1) 1 tdata(end) 0.5 10]);
                ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
                [result,gof] = fit(time(1:end)', tdata(1:end),ft);
                T2f(i,j,k)=min(result.b,result.e);
                T2s(i,j,k)=max(result.b,result.e);
                R2bi(i,j,k)=gof.rsquare;     % R^2 of the biexponential fit

                
                fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[Inf,100, Inf],...
                   'StartPoint',[tdata(1) 5 tdata(end)]);
                ft = fittype('c+a*exp(-x/b)','options',fo); 
                [result,gof] = fit(time(1:end)', tdata(1:end),ft);
                T2(i,j,k)=result.b;
                R2mono(i,j,k)=gof.rsquare;   % R^2 of the mono-exponential fit
            else
                T2f(i,j,k)=0;
                T2s(i,j,k)=0;
                T2(i,j,k)=0;
                R2bi(i,j,k)=0;
                R2mono(i,j,k)=0;
            end

       end
    end
end

%% write the per-voxel maps as AFNI datasets (same cell as T2estimate_scan402_418.m; needed by alignT2toDTI627.m)
[anat,info]=BrikLoad('allave+orig.BRIK');   % grid of the 142-160 series (16x16x16 @ 0.4 mm)
delete(['T2*orig*']);
opt.prefix=['T2f'];
WriteBrik(T2f,info,opt);
opt.prefix=['T2s'];
WriteBrik(T2s,info,opt);
opt.prefix=['T2'];
WriteBrik(T2,info,opt);
opt.prefix=['T2R2bi'];      % R^2 of the biexponential fit per voxel (gof.rsquare); used for the R2 rows of the supp tables
WriteBrik(R2bi,info,opt);
opt.prefix=['T2R2mono'];    % R^2 of the mono-exponential fit per voxel
WriteBrik(R2mono,info,opt);

%%
anat=squeeze(BrikLoad('underlay_164+orig.BRIK'));

% maskt=BrikLoad('COPY_20250101Na.13.an200epi.resamp+orig');
% maskt(:,80:end,:)=0;
% mask=BrikLoad('ROI+orig.BRIK');

mask=BrikLoad('ROI+orig.BRIK');
% mask(mask<12000)=0;
% mask(mask>0)=1;

maskt=mask;
sl=7:8;
% for sl=10%:size(data,3)
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
colormap('jet')
% clim([0 0.5]);
%
alphadata =1*min(maskt(:,:,sl),[],3)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2*short (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2short_map',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2short_map',num2str(sl),'.eps'],'-depsc','-r300');


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
colormap('jet')
% clim([0 55]);
% alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2*long (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2long_map',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2long_map',num2str(sl),'.eps'],'-depsc','-r300');

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
colormap('jet')
% clim([0 55]);
% alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2* (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2estimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2estimate',num2str(sl),'.eps'],'-depsc','-r300');
% end

% %%
% figure;
% 
% 
% 
%  fo = fitoptions('Method','NonlinearLeastSquares',...
%                    'Lower',[0,0,0,0,5],...
%                    'Upper',[Inf,5, Inf,Inf,100],...
%                    'StartPoint',[tdata(1) 0.1 tdata(end) tdata(1) 5]);
% ft = fittype('c+a*exp(-x/b)+d*exp(-x/e)','options',fo); 
% result1 = fit(repmat(time(1:end)',[size(Tdata1,2),1]), Tdata1(:),ft);
% 
% result2 = fit(repmat(time(1:end)',[size(Tdata2,2),1]), Tdata2(:),ft);
% 
% plot(result1)
% hold on
% plot(result2,'b')
% scatter(time,Tdata1,'r','filled');
% 
% scatter(time,Tdata2,'b','filled');
% title({['Cortex T2_f=',num2str(result1.b),'ms T2_s=',num2str(result1.e),'ms'] ...
%     ['CSF T2_f=',num2str(result2.b),'ms T2_s=',num2str(result2.e),'ms'] })
% fontsize(gcf,12,"points")
% ylabel('SNR')
% xlabel('TE (ms)')
% legend({'Cortex fitted result' 'CSF fitted result'});
% print(gcf,[pwd,'/allT2estimate.jpg'],'-djpeg','-r300');
% print(gcf,[pwd,'/allT2estimate.eps'],'-depsc','-r300');
% eee
% 
% %%
% sl=9;
% 
% figure;imagesc(T2(:,:,sl)')
% tdata=squeeze(data(9,10,sl,:));
% fo = fitoptions('Method','NonlinearLeastSquares',...
%    'Lower',[0,0,0],...
%    'Upper',[Inf,50, Inf],...
%    'StartPoint',[tdata(1) 5 tdata(end)]);
% ft = fittype('c+a*exp(-x/b)','options',fo); 
% result = fit(time(1:end)', tdata(1:end),ft);
% result.b
% %
% figure(Position=[100 100 400 300]);
% plot(result,time(1:end)', tdata(1:end))
% fontsize(gcf,12,"points")
% ylabel('SNR')
% xlabel('TE (ms)')
% title(['T2=',num2str(result.b),'ms'])