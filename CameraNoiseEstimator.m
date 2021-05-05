clear all;close all;clc

ROI_number = 10000; %10000 300
subname ='.\roi_';%%row_roi_ %%col_roi_ %% roi_
path1 ='.\Fixde_ROI_QERROR _EachPixel\Save\Data\A4C920PointGrey\cam0';
roi_r  = csvread([path1 subname 'r.csv']);
roi_g  = csvread([path1 subname 'g.csv']);
roi_b  = csvread([path1 subname 'b.csv']);


VarR = getMartrVar(roi_r,1);
VarG = getMartrVar(roi_g,1);
VarB = getMartrVar(roi_b,1);

MeanR = getMartrMean(roi_r,1);
MeanG = getMartrMean(roi_g,1);
MeanB = getMartrMean(roi_b,1);

PixelRange = 1:ROI_number;
figure('Renderer', 'painters', 'Position', [10 10 1200 900])
subplot(3,3,[1 2])
yyaxis right
plot(PixelRange,MeanR(1:ROI_number),'-k');ylabel('mean');title(['Var mean:' num2str(mean(VarR(1:ROI_number)))]);
yyaxis left
plot(PixelRange,VarR(1:ROI_number),'-r');ylabel('var');xlabel('pixel index');%xlim([50 60])
subplot(3,3,[4 5])
yyaxis right
plot(PixelRange,MeanG(1:ROI_number),'-k');ylabel('mean');title(['Var mean:' num2str(mean(VarG(1:ROI_number)))]);
yyaxis left
plot(PixelRange,VarG(1:ROI_number),'-g');ylabel('var');xlabel('pixel index');%xlim([50 60])
subplot(3,3,[7 8])
yyaxis right
plot(PixelRange,MeanB(1:ROI_number),'-k');ylabel('mean');title(['Var mean:' num2str(mean(VarB(1:ROI_number)))]);
yyaxis left
plot(PixelRange,VarB(1:ROI_number),'-b');ylabel('var');xlabel('pixel index');%xlim([50 60])
subplot(3,3,[3 6 9])
VarDataSet =[VarR,VarG,VarB];
boxplot(VarDataSet,'labels',{'VarR','VarG','VarB'});
saveas(gcf,[subname 'PixelMeanVar.bmp'])

% figure()
% plot(abs(fft(VarR-mean(VarR),512)));
% figure()
% plot(abs(fft(VarR-mean(VarR),512)));
% figure()
% plot(abs(fft(MeanB-mean(MeanB),512)));
% figure()
% plot(abs(fft(MeanB-mean(MeanB),512)));
function VarData = getMartrVar(Matric,flag)
[m,n]=size(Matric);
if( flag)
    VarData = zeros(n,1);
    for idx = 1:n
        VarData(idx)=var(Matric(:,idx));
    end
else    
    VarData = zeros(m,1);
    for idx = 1:m
        VarData(idx)=var(Matric(idx,:));
    end
end
end
function MeanData = getMartrMean(Matric,flag)
[m,n]=size(Matric);
if( flag)
    MeanData = zeros(n,1);
    for idx = 1:n
        MeanData(idx)=mean(Matric(:,idx));
    end
else    
    MeanData = zeros(m,1);
    for idx = 1:m
        MeanData(idx)=mean(Matric(idx,:));
    end
end
end