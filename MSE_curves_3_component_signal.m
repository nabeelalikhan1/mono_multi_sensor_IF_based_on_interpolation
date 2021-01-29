clc
clear all;
close all
SampFreq = 256;
addpath('D:\tfsa_5-5\windows\win64_bin');
n=0:255;
N=256;
Sig=1*cos(0.00000001*pi*n.^4+0.2*pi*n)+1*cos(-0.00000001*pi*n.^4+0.9*pi*n)+1*cos(-0.00000001*pi*n.^4+0.7*pi*n)+0*exp(1i*n.^4/(70*N^2)+1i*pi*n/(8))+0*exp(1i*n.^4/(81*N^2)-1i*pi*n.^2/N);
IF_O(1,:)=(0.00000001*n.^3*pi*4+0.2*pi)/(2*pi);
IF_O(2,:)=(-0.00000001*n.^3*pi*4+0.9*pi)/(2*pi);
IF_O(3,:)=(-0.00000001*n.^3*pi*4+0.7*pi)/(2*pi);

SigO =Sig;

num=3;
NS=100;
IF_O=2*IF_O.';%/length(IF_O);
% HADTFD BASED
win_length=125;
%for snr=-10:2:10
iiii=0;
delta=4;
L=32*2;
FFT_length=length(Sig);
for snr=0:5:30
    iiii=iiii+1;
    
    for k1=1:NS
        
        Sig=awgn(SigO,snr,'measured');
        
        for kkkkk=0:3
            % ORIGINAL
            
            if kkkkk==0
                [findex] = FAST_IF(Sig,win_length, num, delta,L*1,0,0)*2*SampFreq;
                
            elseif kkkkk==1
                               findex =FASTEST_IF(Sig,win_length, num, delta,L/2,0,0,32,FFT_length)*2*SampFreq;
              %                                                 findex =FASTEST_IF_BSEARH(Sig,win_length, num, delta,L/2,0,0,32,FFT_length)*2*SampFreq;
 
        elseif kkkkk==2
               findex =FASTEST_IF(Sig,win_length, num, delta,L/2,0,0,16,FFT_length)*2*SampFreq;
   %             findex =FASTEST_IF_BSEARH(Sig,win_length, num, delta,L/2,0,0,16,FFT_length)*2*SampFreq;
            else
                                findex =FASTEST_IF(Sig,win_length, num, delta,L/2,0,0,8,FFT_length)*2*SampFreq;
  %findex =FASTEST_IF_BSEARH(Sig,win_length, num, delta,L/2,0,0,8,FFT_length)*2*SampFreq;
            end
            
            
            
            
            msee=0.1*ones(1,num);
            IF=zeros(1,length(Sig));
            dis=0;
            clear c;
            
            for ii=1:num
                
                t=1:SampFreq;
                IF=findex(ii,:)/length(Sig);
                t=t(5:end-5);
                for i=1:num
                    c(i)=sum(abs(IF(t)'-IF_O(t,i)).^2);
                end
                [a1 b1]=min(c);
                if msee(b1)>=a1(1)/length(Sig)
                    msee(b1)=a1(1)/length(Sig);
                end
                if dis==1
                    figure;
                    plot(t,IF(t),'-',t,IF_O(t,b1),'d');
                end
            end 
                if kkkkk==0
                    mse_FAST_IF_1(k1)=mean(msee);
                elseif kkkkk==1
                    mse_FASTER_IF_32_1(k1)=mean(msee);
                elseif kkkkk==2
                    mse_FASTER_IF_16_1(k1)=mean(msee);
                    
                elseif kkkkk==3
                    mse_FASTER_IF_8_1(k1)=mean(msee);
                    
                end
            end
            
        end
        
        
        
        mse_FAST_IF(iiii)=mean(mse_FAST_IF_1);
        mse_FASTER_IF_64(iiii)=mean(mse_FASTER_IF_32_1);
        mse_FASTER_IF_32(iiii)=mean(mse_FASTER_IF_16_1);
        mse_FASTER_IF_16(iiii)=mean(mse_FASTER_IF_8_1);
        
    end
    snr=0:5:30;
    plot(snr, 10*(log10(mse_FAST_IF)),'-rh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_FASTER_IF_64)),'-bh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_FASTER_IF_32)),'-.k+','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_FASTER_IF_16)),'-.y+','linewidth',4);
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
    legend('Step=1','Step=16','Step=8','Step=4');
    
    % hold on;
    % plot(snr, 10*(log10(mse_proposed)),'-.k+','linewidth',4);
    %
    % hold on;
    % plot(snr, 10*(log10(mse_mb)),'-.y+','linewidth',4);
    %
    % hold on;
    % plot(snr, 10*(log10(mse_spec)),'-.g+','linewidth',4);
    %
    
    
    
    %xlabel('Signal to Noise Ratio');
    %ylabel('Mean Square Error (dB)');
    %legend('FAST_IF','Proposed_IF_16','Proposed_IF_8','Proposed_IF_4');
    % axis([min(snr) max(snr)  -50  0])
