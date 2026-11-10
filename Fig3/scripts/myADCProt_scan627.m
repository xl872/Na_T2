
cd(fullfile(fileparts(mfilename('fullpath')),'..','data','DWI_scan627'));   % release layout: run in the dataset folder
clc;clear all;close all;
mask=BrikLoad('ROI+orig.BRIK');
raw='../bruker/';
date='20250801_201306_Mouse_23Naplus1_08012025_1_81';
data=[];
%%
figure
imagesc(mask(:,:,41)); axis xy  % organoid centre ~ data(31,30,41)
%%

for run=[627]
    % data=squeeze(BrikLoad([date,'.',num2str(run),'.an200epi+orig.BRIK']));
    
    [data,info]=BrikLoad(['allave_t+orig.BRIK']);%use averaged data
    % loads a series of data with a number to string conversion
    
    %filter size of 3 (2*ceil(2*sigma)+1
    for i=1:size(data,4)
        data(:,:,:,i)=imgaussfilt3(data(:,:,:,i),0.5);%-mean(bg,'all');
    end
    %soothes out the voxels using 3D gaussian filter with 0.5 sigma for the

    filename = ([raw,date,'/',num2str(run),'/method']);
    fid = fopen(filename,'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    lines=lines{:};
    lineidx=find(strncmp(lines,'##$PVM_DwEffBval',16));
    lineidx1=find(strncmp(lines,'##$PVM_DwGradVec=',16));
    B0=split([lines{lineidx+1:lineidx1-1}],' ');
    B=cellfun(@str2num,B0);


    sv=size(data);
    V1=reshape(data,[],sv(4));
    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0 -Inf],...
                   'Upper',[0.01 Inf],...
                   'StartPoint',[0.001 B(1)*3]);
    ft = fittype('b*exp(-x*a)','options',fo);

    ADC=zeros(sv(1:3));
    for i=1:size(V1,1)
        if mask(i)
        result = fit(B, V1(i,1:end)',ft);
        ADC(i)=result.a;
        end
    end
    ADC=reshape(ADC,sv(1:3));
    opt.prefix=['ADC'];
    delete([opt.prefix,'*']);
    WriteBrik(ADC*1000000,info,opt);
    
    tdata=mean(V1(find(mask),:),1);
    result=fit(B,tdata(1:end)',ft);
    figure(Position=[100 100 400 300]);
    plot(result,B(1:end)', tdata(1:end))
    fontsize(gcf,12,"points")
    ylabel('Signal Intensity')
    xlabel('B (s/mm^2)')
    title(['ADC=',num2str(result.a),'mm^2/s'])
    print(gcf,[pwd,'/AveADCestimate.jpg'],'-djpeg','-r300');
    print(gcf,[pwd,'/AveADCestimate.eps'],'-depsc','-r300');
end
%%
%%
anat=squeeze(BrikLoad('allave+orig.BRIK'));

% maskt=BrikLoad('COPY_20250101Na.13.an200epi.resamp+orig');
% maskt(:,80:end,:)=0;
% mask=BrikLoad('ROI+orig.BRIK');

mask=BrikLoad('ROI+orig.BRIK');
% mask(mask<12000)=0;
% mask(mask>0)=1;

maskt=mask;
sl=41:41;  % 627: coronal (round-tube) slice through the organoid; z = tube axis, 81 slices
% for sl=10%:size(data,3)
figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
B=mean(anat(:,:,sl),3)';
B=repmat(B,1,1,3)/max(B(:));
hB = image(B);%axis xy
axis image off;
hold on;
F=mean(ADC(:,:,sl),3)';
% F=imresize(F,[size(B,1),size(B,2)]);
hF = imagesc(F);%axis xy
colormap('jet')
% clim([0 0.5]);
%
alphadata =1*min(maskt(:,:,sl),[],3)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'ADC (mm^2/s)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/ADCestimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/ADCestimate',num2str(sl),'.eps'],'-depsc','-r300');
