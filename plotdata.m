clear all;close all;clc
data_snr  = csvread('col\Faker_G_snr.csv');
ROI_number = 7380;
ROI_index = 0:ROI_number-1;
plot(ROI_index,data_snr);hold on
% data_snr  = csvread('row\Faker_R_snr.csv');
% ROI_number = 7380;
% ROI_index = 0:ROI_number-1;
% plot(ROI_index,data_snr);hold on

roi_step =50;
roi_range = [1:roi_step:10000]; % SNR is in dB scale
simu_snr = csvread('Avg_snr.csv');
plot(roi_range,simu_snr);
xlabel('pixel number');xlim([0 7380]) 
ylabel('SNR dB');
legend('fake face','Simulation')
set(findall(gcf,'Type','line'),'LineWidth',2)
saveas(gcf,'compare.bmp')