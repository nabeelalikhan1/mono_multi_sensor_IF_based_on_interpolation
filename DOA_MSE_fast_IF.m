close all;
clear all;
N_sensors=2;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('E:\Published Papers\DOA estimation of intersecting components 2018\Matlab code');
%addpath('E:\Published Papers\DOA ESTIMATION VITERBI\Multi-sensor IF estimation code');

%crossing componentsi8

s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.5*n-0.45*n.^3/(128*128*3)));

perc=0.4;

IF_O(1,:)=0.05+0.45*3*n.^2/(128*128*3);
IF_O(2,:)=0.5-0.45*3*n.^2/(128*128*3);
IF_O=IF_O.';
%s = [(s1.') (s2.') (s3.')];
s = [(s1.') (s2.') ];


n_sources=2;
N_C=2;
s_orig=s;

% set mixing matrix A
%theta = [15,30,50]*pi/180;   % sensor separation angles in radians
%        theta = [-12,0,12]*pi/180;   % sensor separation angles in radians
theta = [-10,10]*pi/180;   % sensor separation angles in radians
LL=100;
index=0;
delta=2;
num=n_sources;
for SNR=-10:2:10
    %for SNR=-2:2:0
    for ii=1:LL
        
        
        A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
        
        
        X = A*s.';                             % mixed source
        theta9=round(theta *180/pi);
        % generate noise
        
        sigma = 10^(-SNR/20);
        w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
        
        X=X+w;
        
        
        %tic
        for k=0:3
        
            switch k 
                case 0
                    [IFF,ss] = Multi_Sensor_FASTEST_IF(X,N_sensors,65, n_sources, delta,32,0,0,32,length(X));
                    
                case 1
                    
                    [IFF,ss] = Multi_Sensor_FASTEST_IF(X,N_sensors,65, n_sources, delta,32,0,0,16,length(X));
                    
                case 2
                    
                    [IFF,ss] = Multi_Sensor_FASTEST_IF(X,N_sensors,65, n_sources, delta,32,0,0,8,length(X));
                    
                case 3
                    
                    [IFF,ss] = Multi_Sensor_FAST_IF(X,N_sensors,65, n_sources, delta,64,0,0);
                    
            end
                    
        
        
        %toc
        clear a;
        for iii=1:n_sources
            for jjj=1:N_sensors
                a(jjj,:)=ss(jjj,iii,:);
            end
            theta1=-90:0.5:90;
            
            p=TMMUSIC(cov(a.'), 2, N_sensors, 1, 1, theta1');
           % figure; plot(p);
            [x,y]=max(p);
            y1(iii)=y(1)/2;
        end
        
        y1=y1-90;
        
        
                    msee=0.1*ones(1,num);

        for ii22=1:num
            
            t=1:128;
            IF=IFF(ii22,:);%/length(X);
            t=t(5:end-5);
            for i=1:num
                c(i)=sum(abs(IF(t)'-IF_O(t,i)).^2);
            end
            [a1, b1]=min(c);
            if msee(b1)>=a1(1)/length(X)
                msee(b1)=a1(1)/length(X);
            end
            
        end
        
        switch k
            case 0
                mseeIF_step_size_32(ii)=mean(msee);
                mmssee_step_size_32(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);
            case 1
                mseeIF_step_size_16(ii)=mean(msee);
                mmssee_step_size_16(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);
            case 2
                mseeIF_step_size_8(ii)=mean(msee);
                mmssee_step_size_8(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);
            case 3
                mseeIF_step_size_1(ii)=mean(msee);
                mmssee_step_size_1(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);
            
        end
        
        end
    end
    index=index+1;
    %mean(mmssee)
    snr_mse_step_size_8(index)=mean(mmssee_step_size_8);
    snr_mse_step_size_32(index)=mean(mmssee_step_size_32);
    snr_mse_step_size_16(index)=mean(mmssee_step_size_16);
    snr_mse_step_size_1(index)=mean(mmssee_step_size_1);
    
    IF_mse_step_size_8(index)=mean(mseeIF_step_size_8);
    IF_mse_step_size_32(index)=mean(mseeIF_step_size_32);
    IF_mse_step_size_16(index)=mean(mseeIF_step_size_16);
    IF_mse_step_size_1(index)=mean(mseeIF_step_size_1);
    
    
end

%SNR=-10:2:0;
SNR=-10:2:10;
% FROM SIMULATION DONE in "Novel direction of arrival estimation using
% spatial adaptive"
%snr_mse_sadtfd=[0.4418    0.0994    0.0615    0.0414    0.0272    0.0204];
plot(SNR,10*(log10(snr_mse_step_size_1)),'--md','linewidth',2);
hold on;
plot(SNR,10*(log10(snr_mse_step_size_8)),'r','linewidth',3);
hold on;
plot(SNR,10*(log10(snr_mse_step_size_16)),'g','linewidth',3);
hold on;
%plot(SNR,10*(log10(snr_mse_post_proc)),'y','linewidth',2);
%hold on;
plot(SNR,10*(log10(snr_mse_step_size_32)),'b:','linewidth',3);

xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('Step size 1','Step size 8','Step size 16','Step size 32');
title('DOA estimation accuracy');
%legend('The Proposed Method','Time-frequency Music','DOA based on IF estimation using ridge tracking');


SNR=-10:2:10;
figure;
% FROM SIMULATION DONE in "Novel direction of arrival estimation using
% spatial adaptive"
%snr_mse_sadtfd=[0.4418    0.0994    0.0615    0.0414    0.0272    0.0204];
plot(SNR,10*(log10(IF_mse_step_size_1)),'--md','linewidth',2);
hold on;
plot(SNR,10*(log10(IF_mse_step_size_8)),'r','linewidth',3);
hold on;
plot(SNR,10*(log10(IF_mse_step_size_16)),'g','linewidth',3);
hold on;
%plot(SNR,10*(log10(snr_mse_post_proc)),'y','linewidth',2);
%hold on;
plot(SNR,10*(log10(IF_mse_step_size_32)),'b:','linewidth',3);

xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('Step size 1','Step size 8','Step size 16','Step size 32');
title('IF estimation accuracy');


