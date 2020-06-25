% ********************* BPSK modulation and de-modulation ****************%
clc;
clear all;
close all;
% ********************* Define transmitted signal *************************
N=10; % Number of bits , size of transmitted signal x_inp=[x_1 x_2 ... x_N] 
x_inp= [1 0 1 1 0 1 1 1 0 1];  % binary signal 0 or 1 % message to be transmitted                               
Tb=0.0001; % bit period (second)   
% ********************* Represent input signal as digital signal ****
x_bit=[]; 
nb=100; % bbit/bit
for n=1:1:N   % 
    if x_inp(n)==1;  % 
       x_bitt=ones(1,nb);
    else x_inp(n)==0;
        x_bitt=-x_bitt;
    end
     x_bit=[x_bit x_bitt];
end
t1=Tb/nb:Tb/nb:nb*N*(Tb/nb); % time of the signal 
f1 = figure(1);
set(f1,'color',[1 1 1]);
subplot(4,1,1);
plot(t1,x_bit,'lineWidth',2);grid on;
axis([ 0 Tb*N -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Input signal');
% % ********************* Define BFSK Modulation ****************************
% Ac=5;  % Amplitude of carrier signal
% mc=4;  % fc>>fs fc=mc*fs fs=1/Tb
% fc=mc*(1/Tb); % carrier frequency for bit 1
% fi1=0; % carrier phase for bit 1
% fi2=pi; % carrier phase for bit 0
% t2=Tb/nb:Tb/nb:Tb;                 
% t2L=length(t2);
% x_mod=[];
% for (i=1:1:N)
%     if (x_inp(i)==1)
%         x_mod0=Ac*cos(2*pi*fc*t2+fi1);%modulation signal with carrier signal 1
%     else
%         x_mod0=Ac*cos(2*pi*fc*t2+fi2);%modulation signal with carrier signal 2
%     end
%     x_mod=[x_mod x_mod0];
% end
% t3=Tb/nb:Tb/nb:Tb*N;
% subplot(3,1,2);
% plot(t3,x_mod);
% xlabel('Time(sec)');
% ylabel('Amplitude(volt)');
% title('BPSK modulated signal ');
% ********************* Transmitted signal x ******************************
% x=x_mod;
% % ********************* Channel model h and w *****************************
% h=1; % Fading 
% w=0; % Noise
% % ********************* Received signal y *********************************
% y=h.*x+w;
% Fs=nb/Tb;
% figure
% [Pxx ff]=pwelch(y);
% plot(ff/pi*Fs/2/1e3,db(abs(Pxx))),grid
% ylabel('Power');
% xlabel(' Frequency(kHz)');
% title(' PSD Çýktýsý ');
% 
% nfft=4096;
% X = fft(y,nfft);
% X = X(1:nfft/2);
% % Take the magnitude of fft of x
% signalFFT = abs(X);
% f = (0:nfft/2-1)*Fs/nfft;
% figure(1),subplot(5,1,1);plot(f/1e3,signalFFT,'b');grid
% xlabel('Frequency(kHz)');
% ylabel('Amplitude');
% title('PSK iþareti Frekans Spektrumu')
% figure(1),subplot(5,1,2);plot(f/1e3,signalFFT,'b');grid
% xlabel('Frequency(kHz)');
% ylabel('Amplitude');
% title('PSK iþareti Frekans Spektrumu')

data=[1 0 1 1 0 1 1 1 0 1]; % information
%Number_of_bit=1024;
%data=randint(Number_of_bit,1);
% figure(1)
% stem(data, 'linewidth',3), grid on;
% title('  Information before Transmiting ');
% axis([ 0 11 0 1.5]);
data_NZR=2*data-1; % Data Represented at NZR form for QPSK modulation
s_p_data=reshape(data_NZR,2,length(data)/2);  % S/P convertion of data
br=10.^4; %Let us transmission bit rate  1000000
f=4*br; % minimum carrier frequency
T=1/br; % bit duration
t=T/100:T/100:T; % Time vector for one bit information
% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulatio  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
Tb=0.0001;
N=10;
y=[];
y_in=[];
y_qd=[];
for(i=1:length(data)/2)
    y1=s_p_data(1,i)*cos(2*pi*f*t)*0.707; % inphase component
    y2=s_p_data(2,i)*sin(2*pi*f*t)*0.707 ;% Quadrature component 90 derece
    y_in=[y_in y1]; % inphase signal vector
    y_qd=[y_qd y2]; %quadrature signal vector
    y=[y y1+y2]; % modulated signal vector
end
Tx_sig=y; % transmitting signal after modulation
tt=T/100:T/100:(T*length(data))/2;
% figure(2)
subplot(4,1,2);
plot(tt,y_in,'linewidth',1), grid on;
axis([ 0 Tb*N -1.5 1.5]);
title(' I Kanalý');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

subplot(4,1,3);

plot(tt,y_qd,'linewidth',1), grid on;
axis([ 0 Tb*N -1.5 1.5]);
title(' Q Kanalý');
xlabel('time(sec)');
ylabel(' amplitude(volt0');
subplot(4,1,4);

plot(tt,Tx_sig,'r','linewidth',1), grid on;
title('I+Q , QPSK Modulason Sonucu');
axis([ 0 Tb*N -1.5 1.5]);
xlabel('time(sec)');
ylabel(' amplitude(volt0');


Fs=100/T;
figure;
subplot(3,1,1);
[Pxx ff]=pwelch(Tx_sig);
plot(ff/pi*Fs/2/1e3,db(abs(Pxx))),grid
ylabel('Power');
xlabel(' Frequency(kHz)');
title(' PSD Çýktýsý ');

nfft=4096;
X = fft(Tx_sig,nfft);
X = X(1:nfft/2);
% Take the magnitude of fft of x
signalFFT = abs(X);
f = (0:nfft/2-1)*Fs/nfft;
figure;
subplot(3,1,2);
plot(f/1e3,signalFFT,'b');grid
xlabel('Frequency(kHz)');
ylabel('Amplitude');
title('PSK iþareti Frekans Spektrumu')
