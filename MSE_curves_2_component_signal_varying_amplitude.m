clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;
%Sig1 = 1*exp(1i*(2*pi*(50*t.^1))+1i*(2*pi*(0*t))); %300t»òÕß150t
%Sig2 = 1*exp(1i*(2*pi*(50*t.^2))+1i*(2*pi*(10*t))); %300t»òÕß150t

%Sig3 = exp(1i*(1*pi*(10*t +35*t.^3)));
%Sig4 =1*exp(1i*(1*pi*(110*t -35*t.^3)));
%SigO =1*Sig1 +1*Sig4 +Sig3;
%IF_O(:,1)=100*t/1;
%IF_O(:,2)=100*t/1+10;
%IF_O(:,1)=105*t.^2/2+5;
%IF_O(:,2)=-105*t.^2/2+55;
%IF_O(:,3)=50;


Sig1 = 1*exp(1i*(1*pi*(30*t.^3))+1i*(2*pi*(0*t))); %300t»òÕß150t
Sig2 = 1*exp(1i*(-1*pi*(30*t.^3))+1i*(1*pi*(90*t))); %300t»òÕß150t

Sig3 = exp(1i*(1*pi*(20*t +30*t.^3)));
Sig =1*Sig1 +0*Sig3+0.5*Sig2;
SigO =Sig;
IF_O(:,1)=90*t.^2/2;
IF_O(:,2)=-90*t.^2/2+90/2;
%IF_O(:,3)=90*t.^2/2+10;


%Sig=Sig.*([1:128 128:-1:1]);
num=2;
NS=250;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED
win_length=65;
%for snr=-10:2:10
iiii=0;
delta=2;
L=32*1;
FFT_length=length(Sig);
for snr=0:2:10
    iiii=iiii+1;
    
    for k1=1:NS
        
        Sig=awgn(SigO,snr,'measured');
        
        for kkkkk=0:3
            
            % ORIGINAL
            
            if kkkkk==0
                [findex] = FAST_IF(Sig,win_length, num, delta,L,0,0)*2*SampFreq;
                
            elseif kkkkk==1
                findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,16,length(Sig))*2*SampFreq;
%
      %                     findex =FASTEST_IF_BSEARH(Sig,win_length, num, delta,L,0,0,16,FFT_length)*2*SampFreq;

            elseif kkkkk==2
                findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,8,length(Sig))*2*SampFreq;
         %                                  findex =FASTEST_IF_BSEARH(Sig,win_length, num, delta,L,0,0,8,FFT_length)*2*SampFreq;

            else
            
                findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,4,length(Sig))*2*SampFreq;
                %                                           findex =FASTEST_IF_BSEARH(Sig,win_length, num, delta,L,0,0,4,FFT_length)*2*SampFreq;

                
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
        mse_FASTER_IF_32(iiii)=mean(mse_FASTER_IF_32_1);
        mse_FASTER_IF_16(iiii)=mean(mse_FASTER_IF_16_1);
        mse_FASTER_IF_8(iiii)=mean(mse_FASTER_IF_8_1);
        
    end
    snr=0:2:10;
    plot(snr, 10*(log10(mse_FAST_IF)),'-rh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_FASTER_IF_32)),'-bh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_FASTER_IF_16)),'-.k+','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_FASTER_IF_8)),'-.y+','linewidth',4);
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
    legend('Step size=1','Step size=16','Step size=8','Step size=4');
    
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
