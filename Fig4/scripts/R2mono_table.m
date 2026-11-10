% R2 of the voxel-wise mono-exponential T2* fit (same data prep and fit as
% T2estimate_M1T2.m and T2estimate_M3T2.m), written to R2mono+orig and
% summarized as mean/std within ROI+orig ("Organoid") and ROI_metrics+orig ("Metrics").
cd(fullfile(fileparts(mfilename('fullpath')),'..','..','data'));   % release layout: run in the dataset folder
addpath(genpath(fullfile(getenv('HOME'),'MATLAB Add-Ons','Toolboxes','afni')));
base=pwd;
S(1).dir='M1'; S(1).name='#3'; S(1).pre='2025051802Na.'; S(1).runs=155:173;
S(1).time=[0.88 1 1.2 1.5 1.8 2.4 3 4 5 6 7 8.5 10 15 20 30 40 60 80];
S(2).dir='M3'; S(2).name='#4'; S(2).pre='20250518_144207_Mouse_23Naplus1_05182025_02_1_66.'; S(2).runs=61:79;
S(2).time=[0.68 0.84 1 1.25 1.5 2 3 4 5 6.5 8 10 15 20 30 40 50 65 80];
for s=1:2
    cd(fullfile(base,S(s).dir));
    time=S(s).time; data=[]; BG=[];
    for run=S(s).runs
        [V,info]=BrikLoad([S(s).pre,num2str(run),'.an200epi+orig.BRIK']);
        data0=squeeze(V);
        if s==1, bg=data0(:,[1:4,17:20],[1:4,17:20]); else, bg=data0([1:3,14:16],[1:3,14:16],:); end
        BG=cat(2,BG,bg(:));
        data=cat(4,data,imgaussfilt3(data0,0.5));
    end
    mask=BrikLoad('formask+orig.BRIK'); mask=mask>=10000;
    data_std=std(BG);
    T2=zeros(size(mask)); R2=zeros(size(mask));
    for idx=find(mask)'
        [i,j,k]=ind2sub(size(mask),idx);
        tdata=squeeze(data(i,j,k,:))./data_std';
        fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0],'Upper',[Inf,100, Inf],'StartPoint',[tdata(1) 5 tdata(end)]);
        ft = fittype('c+a*exp(-x/b)','options',fo);
        [result,gof] = fit(time(1:end)', tdata(1:end),ft);
        T2(idx)=result.b; R2(idx)=gof.rsquare;
    end
    T2old=BrikLoad('T2+orig.BRIK');
    fprintf('%s: max |T2 refit - T2+orig| = %.2e ms\n',S(s).name,max(abs(T2(:)-T2old(:))));
    info.BRICK_TYPES=3; opt.OverWrite='y'; opt.prefix='R2mono';
    WriteBrik(R2,info,opt);
    for r={'ROI','ROI_metrics'}
        roi=BrikLoad([r{1},'+orig.BRIK']); v=roi>0;
        fprintf('%s %-11s n=%4d  T2 mono = %.2f +/- %.2f   R2 mono = %.3f +/- %.3f\n',S(s).name,r{1},nnz(v),mean(T2old(v)),std(T2old(v)),mean(R2(v)),std(R2(v)));
    end
end
cd(base);
