clear all;close all;clc
rng('default');

MonteCarloPara.A = 0.5;
MonteCarloPara.N = 512;
MonteCarloPara.T = 1/30; % sec
MonteCarloPara.fo = 1;   % Hz
MonteCarloPara.DC = 0; % DC level
MonteCarloPara.RollingDelay = 0; % delay time

camera_noise = [1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0];
Nrun = 500; % Montecarlo runs
Roi_Step = 50;
Pixel_Range = [50:Roi_Step:10000];
Pixel_Range = [1,Pixel_Range];
% SavePath ='.\SavResultUniformNoise\';
SavePath ='.\SavResultGaussionNoise\';
mkdir(SavePath);

for id_camera_noise =1:length(camera_noise)
    Avg_snr_1 = zeros(length(Pixel_Range),1);
    Avg_snr_2 = zeros(length(Pixel_Range),1);
    Avg_snr_3 = zeros(length(Pixel_Range),1);
    Avg_snr_4 = zeros(length(Pixel_Range),1);       
    Camera_noise = camera_noise(id_camera_noise); 
    parfor id_roi = 1:length(Pixel_Range)           
        Pixel_num = Pixel_Range(id_roi);       
        disp(['Camera noise : ' num2str(Camera_noise) '  |  pixel number :' num2str(id_roi) '/ ' num2str(length(Pixel_Range))]);
        % MonteCarlo
        [Sig_snr_1,Sig_snr_2,Sig_snr_3,Sig_snr_4] = QErrorMonteCarlo(Nrun,Pixel_num,Camera_noise,MonteCarloPara);
        
        Avg_snr_1(id_roi)=mean(Sig_snr_1);
        Avg_snr_2(id_roi)=mean(Sig_snr_2);        
        Avg_snr_3(id_roi)=mean(Sig_snr_3);
        Avg_snr_4(id_roi)=mean(Sig_snr_4);
    end
    
    %plot
    figure(1);
    subplot(211)
    plot(Pixel_Range,Avg_snr_1,'-k');hold on
    plot(Pixel_Range,Avg_snr_2,'-b');
    xlabel('pixel number');title('Normal model');legend('model1','model2');
    ylabel('SNR');
    subplot(212)
    plot(Pixel_Range,Avg_snr_3,'-k');hold on
    plot(Pixel_Range,Avg_snr_4,'-b');
    xlabel('pixel number');title(['Delay model:' num2str(MonteCarloPara.RollingDelay) 'ms']);
    ylabel('SNR');legend('model1','model2');
    
    figure(2);
    subplot(211)
    plot(Pixel_Range,Avg_snr_1,'-k');hold on
    plot(Pixel_Range,Avg_snr_3,'-b');
    xlabel('pixel number');title('model1');legend('Normal','Delay');
    ylabel('SNR');
    subplot(212)
    plot(Pixel_Range,Avg_snr_2,'-k');hold on
    plot(Pixel_Range,Avg_snr_4,'-b');
    xlabel('pixel number');title('model2');legend('Normal','Delay');
    ylabel('SNR');
    
    %save
    csvwrite([SavePath 'AvgSnr_Nomal_Model1_' num2str(camera_noise(id_camera_noise)) '.csv'],Avg_snr_1);
    csvwrite([SavePath 'AvgSnr_Nomal_Model2_' num2str(camera_noise(id_camera_noise)) '.csv'],Avg_snr_2);
    csvwrite([SavePath 'AvgSnr_Delay_Model1_' num2str(camera_noise(id_camera_noise)) '.csv'],Avg_snr_3);
    csvwrite([SavePath 'AvgSnr_Delay_Model2_' num2str(camera_noise(id_camera_noise)) '.csv'],Avg_snr_4);
    saveas(figure(1),[SavePath '1_simulabtion_snr_' num2str(camera_noise(id_camera_noise)) '.bmp'])
    saveas(figure(2),[SavePath '2_simulabtion_snr_' num2str(camera_noise(id_camera_noise)) '.bmp'])
    close all
end


function [Sig_snr_1,Sig_snr_2,Sig_snr_3,Sig_snr_4]= QErrorMonteCarlo(Nrun,Pixel_num,Camera_noise,MonteCarloPara)
    A = MonteCarloPara.A;
    N = MonteCarloPara.N;
    T = MonteCarloPara.T; % sec
    fo = MonteCarloPara.fo; % Hz
    DC = MonteCarloPara.DC; % DC level
    RollingDelay = MonteCarloPara.RollingDelay;
    n = [1:N];    
    Sig_snr_1 = zeros(Nrun,1);  
    Sig_snr_2 = zeros(Nrun,1);   
    Sig_snr_3 = zeros(Nrun,1);  
    Sig_snr_4 = zeros(Nrun,1);   
    for run = 1:Nrun
        x1 = zeros(Pixel_num,N); % model 1      
        x2 = zeros(Pixel_num,N); % model 2
        x3 = zeros(Pixel_num,N); % model 1 delay     
        x4 = zeros(Pixel_num,N); % model 2 delay        
        for idx =1:Pixel_num
         % normal
         % w = sqrt(Camera_noise*12)*rand(1,N);
            w = Camera_noise*randn(1,N);
            anal_x = A*cos(2*pi*fo*T*n) + DC ;
            x1(idx,:) = round(anal_x+ w);
            x2(idx,:) = anal_x+ w;
            % delay
            delay = RollingDelay*ceil(idx/100)/0.001; % ms
            anal_x_delay = A*cos(2*pi*fo*T*n+2*pi*fo*delay) + DC ;
            x3(idx,:) = round(anal_x_delay+ w);
            x4(idx,:) = anal_x_delay + w;
        end
        if Pixel_num ==1
            sig1 = x1;
            sig2 = x2;
            sig3 = x3;
            sig4 = x4;
        else
            sig1 = mean(x1);
            sig2 = mean(x2);
            sig3 = mean(x3);
            sig4 = mean(x4);
        end            
        Sig_snr_1(run) =snr(sig1);
        Sig_snr_2(run) =snr(sig2);
        Sig_snr_3(run) =snr(sig3);
        Sig_snr_4(run) =snr(sig4);
    end
end