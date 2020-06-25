clc;
clear all;
close all;
% ********************* Define transmitted signal *************************
N=10; % Number of bits , size of transmitted signal x_inp=[x_1 x_2 ... x_N]
x_inp=[1 0 1 0 1 1 0 0 1 0];  % binary signal 0 or 1 % message to be transmitted
Tb=0.001; % bit period (second)
% ********************* Represent input signal as digital signal ****
x_bit=[];
nb=100; % bbit/bit
for n=1:1:N   %
    if x_inp(n)==1;  %
        x_bitt=ones(1,nb);
    else x_inp(n)==0;
        x_bitt=zeros(1,nb);
    end
    x_bit=[x_bit x_bitt];
end
t1=Tb/nb:Tb/nb:nb*N*(Tb/nb); % time of the signal
% f1 = figure(1);
% set(f1,'color',[1 1 1]);
% subplot(5,1,1);
% plot(t1,x_bit,'lineWidth',2);grid on;
% axis([ 0 Tb*N -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Input signal as digital signal');
% ********************* Define BFSK Modulation ****************************
Ac=2;  % Amplitude of carrier signal
mc1=6;  % fc>>fs fc=mc*fs fs=1/Tb

SG =5; % zero crossing of cross-corelation diagram (also run for 2,3,4, and 5 values)

fc1=mc1*(1/Tb) % carrier frequency for bit 1
fc0=fc1-SG*(1/(2*Tb)) % carrier frequency for bit 0


t2=Tb/nb:Tb/nb:Tb;                 
t2L=length(t2);

c1=Ac*cos(2*pi*fc1.*t2);
c0=Ac*cos(2*pi*fc0.*t2);

% check the signal orthogonality 
zeynep=c1.*c0;    %% s1(t) ve s(0)t
integral=trapz(Tb,zeynep) %Tb suresinde integral
rhoo=(1/Tb)*(integral)% rho boyle hesaplanýyor


% figure(2),plot(t2,c1)
% hold on,plot(t2,c0,'--r'),grid
% legend('fc_1','fc_0')
t2L=length(t2);
x_mod=[];
for i=1:N
    if x_inp(i)==1
        x_mod0=Ac*cos(2*pi*fc1*t2);%modulation signal with carrier signal 1
    else
        x_mod0=Ac*cos(2*pi*fc0*t2);%modulation signal with carrier signal 2
    end
    x_mod=[x_mod x_mod0];
end
t3=Tb/nb:Tb/nb:Tb*N;
figure(1),subplot(5,1,1);
plot(t3,x_mod);
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
title('SG=1 Signal of BFSK modulation ');

% ********************* Transmitted signal x ******************************
x=x_mod;
% ********************* Channel model h and w *****************************
h=1; % Fading
w=0; % Noise
% ********************* Received signal y *********************************
y=h.*x+w;
% ********************* Define BFSK Demodulation **************************
Fs=nb/Tb;
figure
[Pxx ff]=pwelch(y);
% plot(ff/pi*Fs/2/1e3,db(abs(Pxx))),grid
% ylabel('Power');
% xlabel(' Frequency(kHz)');
% title('PSD of FSK Signal');

nfft=4096;
X = fft(y,nfft);
X = X(1:nfft/2);
% Take the magnitude of fft of x
signalFFT = abs(X);
f = (0:nfft/2-1)*Fs/nfft;
% figure(1),subplot(4,1,3);plot(f/1e3,signalFFT,'b');grid
% xlabel('Frequency(kHz)');
% ylabel('Amplitude');
%  Demod
% Part B: Demodulation
bandpass_mod=[]; %BPF yazdýrmak için blgileri tutturacaðýz
y_dem=[]; %iþaretleri ayrý ayrý tutacaðýz
t=Tb/nb:Tb/nb:Tb;
c_dem1=cos(2*pi*fc1*t); % carrier siignal for information 1
c_dem2=cos(2*pi*fc0*t); % carrier siignal for information 0
dogrultucu_mod=[];   %dogrultucu bilgisini yazdýrmak için
rend_mod=[];    %LPF filte çýkýþý için
for n=t2L:t2L:length(y)
  
  if (x_inp(10)==1)
  y_dem1=c_dem2.*y((n-(t2L-1)):n);  
  
  
  dogrultucu1=abs(y_dem1);   %dogrultucu
  [B,A] = butter(6,(fc1+500)/100000,'low');  %%LPF
  rend=filter(B,A,dogrultucu1);    %%LPF
  else
      y_dem1=c_dem1.*y((n-(t2L-1)):n);
      dogrultucu1=abs(y_dem1);
  [B,A] = butter(6,(fc0-500)/100000,'low');
  rend=filter(B,A,dogrultucu1);
  end
      
  
  if(rend<=1.137);      % % logic level = (Ac)/2
    a=0;
  else
    a=1;
  end
  y_dem=[y_dem a];
  rend_mod=[rend_mod rend];
  dogrultucu_mod=[dogrultucu_mod dogrultucu1];
  bandpass_mod=[bandpass_mod y_dem1];
end
x_out=y_dem; % output signal;

% *************** Represent output signal as digital signal ***************
xx_bit=[];
for n=1:length(x_out);
    if x_out(n)==1;
       xx_bitt=ones(1,nb);
    else x_out(n)==0;
        xx_bitt=zeros(1,nb);
    end
     xx_bit=[xx_bit xx_bitt];
end
t4=Tb/nb:Tb/nb:nb*length(x_out)*(Tb/nb);
subplot(4,1,4)
plot(t4,xx_bit,'LineWidth',2);grid on;
axis([ 0 Tb*length(x_out) -0.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('SG=5 Output signal as digital signal');

subplot(4,1,3)
plot(t4,rend_mod,'LineWidth',2);grid on;
axis([ 0 Tb*length(x_out) -0.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('SG=5 LPF çýkýþý');
subplot(4,1,2)
plot(t4,dogrultucu_mod,'LineWidth',2);grid on;
axis([ 0 Tb*length(x_out) -0.5 3]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('SG=5 Doðrultucu çýkýþý');

subplot(4,1,1)
plot(t4,bandpass_mod,'LineWidth',2);grid on;
axis([ 0 Tb*length(x_out) -0.5 3]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('SG=5 BPF çýkýþý');
% **************************** end of program *****************************