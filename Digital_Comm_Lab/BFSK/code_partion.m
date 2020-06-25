% ********************* BFSK modulation and de-modulation ****************%
clc;
clear all;
close all;
% ********************* Define transmitted signal *************************
N=10; % Number of bits , size of transmitted signal x_inp=[x_1 x_2 ... x_N]
x_inp= [1 0 1 0 1 1 0 0 1 0]  % binary signal 0 or 1 % message to be transmitted
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
mc1=5;  % fc>>fs fc=mc*fs fs=1/Tb

SG =1; % zero crossing of cross-corelation diagram (also run for 2,3,4, and 5 values)

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


figure(2),plot(t2,c1)
hold on,plot(t2,c0,'--r'),grid
legend('fc_1','fc_0')
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
title('SG=5 Signal of BFSK modulation ');

% % SG =2; % zero crossing of cross-corelation diagram (also run for 2,3,4, and 5 values)
% % 
% fc1=mc1*(1/Tb) % carrier frequency for bit 1
% fc0=fc1-SG*(1/(2*Tb)) % carrier frequency for bit 0
% 
% 
% t2=Tb/nb:Tb/nb:Tb;                 
% t2L=length(t2);
% 
% c1=Ac*cos(2*pi*fc1.*t2);
% c0=Ac*cos(2*pi*fc0.*t2);
% 
% % check the signal orthogonality 
% zeynep=c1.*c0;    %% s1(t) ve s(0)t
% integral=trapz(Tb,zeynep) %Tb suresinde integral
% rhoo=(1/Tb)*(integral)% rho boyle hesaplanýyor
% 
% 
% % figure(2),plot(t2,c1)
% % hold on,plot(t2,c0,'--r'),grid
% % legend('fc_1','fc_0')
% t2L=length(t2);
% x_mod=[];
% for i=1:N
%     if x_inp(i)==1
%         x_mod0=Ac*cos(2*pi*fc1*t2);%modulation signal with carrier signal 1
%     else
%         x_mod0=Ac*cos(2*pi*fc0*t2);%modulation signal with carrier signal 2
%     end
%     x_mod=[x_mod x_mod0];
% end
% t3=Tb/nb:Tb/nb:Tb*N;
% figure(1),subplot(5,2,3);
% plot(t3,x_mod);
% xlabel('Time(sec)');
% ylabel('Amplitude(volt)');
% title('SG=2 Signal of BFSK modulation ');
% SG =3; % zero crossing of cross-corelation diagram (also run for 2,3,4, and 5 values)
% 
% fc1=mc1*(1/Tb) % carrier frequency for bit 1
% fc0=fc1-SG*(1/(2*Tb)) % carrier frequency for bit 0
% 
% 
% t2=Tb/nb:Tb/nb:Tb;                 
% t2L=length(t2);
% 
% c1=Ac*cos(2*pi*fc1.*t2);
% c0=Ac*cos(2*pi*fc0.*t2);
% 
% % check the signal orthogonality 
% zeynep=c1.*c0;    %% s1(t) ve s(0)t
% integral=trapz(Tb,zeynep) %Tb suresinde integral
% rhoo=(1/Tb)*(integral)% rho boyle hesaplanýyor
% 
% 
% % figure(2),plot(t2,c1)
% % hold on,plot(t2,c0,'--r'),grid
% % legend('fc_1','fc_0')
% t2L=length(t2);
% x_mod=[];
% for i=1:N
%     if x_inp(i)==1
%         x_mod0=Ac*cos(2*pi*fc1*t2);%modulation signal with carrier signal 1
%     else
%         x_mod0=Ac*cos(2*pi*fc0*t2);%modulation signal with carrier signal 2
%     end
%     x_mod=[x_mod x_mod0];
% end
% t3=Tb/nb:Tb/nb:Tb*N;
% figure(1),subplot(5,2,5);
% plot(t3,x_mod);
% xlabel('Time(sec)');
% ylabel('Amplitude(volt)');
% title('SG=3 Signal of BFSK modulation ');
% SG =4; % zero crossing of cross-corelation diagram (also run for 2,3,4, and 5 values)
% 
% fc1=mc1*(1/Tb) % carrier frequency for bit 1
% fc0=fc1-SG*(1/(2*Tb)) % carrier frequency for bit 0
% 
% 
% t2=Tb/nb:Tb/nb:Tb;                 
% t2L=length(t2);
% 
% c1=Ac*cos(2*pi*fc1.*t2);
% c0=Ac*cos(2*pi*fc0.*t2);
% 
% % check the signal orthogonality 
% zeynep=c1.*c0;    %% s1(t) ve s(0)t
% integral=trapz(Tb,zeynep) %Tb suresinde integral
% rhoo=(1/Tb)*(integral)% rho boyle hesaplanýyor
% 
% 
% % figure(2),plot(t2,c1)
% % hold on,plot(t2,c0,'--r'),grid
% % legend('fc_1','fc_0')
% t2L=length(t2);
% x_mod=[];
% for i=1:N
%     if x_inp(i)==1
%         x_mod0=Ac*cos(2*pi*fc1*t2);%modulation signal with carrier signal 1
%     else
%         x_mod0=Ac*cos(2*pi*fc0*t2);%modulation signal with carrier signal 2
%     end
%     x_mod=[x_mod x_mod0];
% end
% t3=Tb/nb:Tb/nb:Tb*N;
% figure(1),subplot(5,2,7);
% plot(t3,x_mod);
% xlabel('Time(sec)');
% ylabel('Amplitude(volt)');
% title('SG=4 Signal of BFSK modulation ');
% SG =5; % zero crossing of cross-corelation diagram (also run for 2,3,4, and 5 values)
% 
% fc1=mc1*(1/Tb) % carrier frequency for bit 1
% fc0=fc1-SG*(1/(2*Tb)) % carrier frequency for bit 0
% 
% 
% t2=Tb/nb:Tb/nb:Tb;                 
% t2L=length(t2);
% 
% c1=Ac*cos(2*pi*fc1.*t2);
% c0=Ac*cos(2*pi*fc0.*t2);
% 
% % check the signal orthogonality 
% zeynep=c1.*c0;    %% s1(t) ve s(0)t
% integral=trapz(Tb,zeynep) %Tb suresinde integral
% rhoo=(1/Tb)*(integral)% rho boyle hesaplanýyor
% 
% 
% % figure(2),plot(t2,c1)
% % hold on,plot(t2,c0,'--r'),grid
% % legend('fc_1','fc_0')
% t2L=length(t2);
% x_mod=[];
% for i=1:N
%     if x_inp(i)==1
%         x_mod0=Ac*cos(2*pi*fc1*t2);%modulation signal with carrier signal 1
%     else
%         x_mod0=Ac*cos(2*pi*fc0*t2);%modulation signal with carrier signal 2
%     end
%     x_mod=[x_mod x_mod0];
% end
% t3=Tb/nb:Tb/nb:Tb*N;
% figure(1),subplot(5,2,9);
% plot(t3,x_mod);
% xlabel('Time(sec)');
% ylabel('Amplitude(volt)');
% title('SG=5 Signal of BFSK modulation ');

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
plot(ff/pi*Fs/2/1e3,db(abs(Pxx))),grid
ylabel('Power');
xlabel(' Frequency(kHz)');
title('SG=2 PSD of FSK Signal');

nfft=4096;
X = fft(y,nfft);
X = X(1:nfft/2);
% Take the magnitude of fft of x
signalFFT = abs(X);
f = (0:nfft/2-1)*Fs/nfft;
figure(1),subplot(5,1,2);plot(f/1e3,signalFFT,'b');grid
xlabel('Frequency(kHz)');
ylabel('Amplitude');
title('Frekans Spektrumu')
% %  Demod
% Part B: Demodulation

% y_dem=[];
% t=Tb/nb:Tb/nb:Tb;
% c_dem1=cos(2*pi*fc1*t); % carrier siignal for information 1
% c_dem2=cos(2*pi*fc0*t); % carrier siignal for information 0
% for n=t2L:t2L:length(y)
%   
%   
%   y_dem1=c_dem1.*y((n-(t2L-1)):n);
%   y_dem2=c_dem2.*y((n-(t2L-1)):n);
%   t4=Tb/nb:Tb/nb:Tb;
%   z1=trapz(t4,y_dem1);  % intregation 
%   z2=trapz(t4,y_dem2);  % intregation 
%   A_dem1=round(2*z1/Tb);
%   A_dem2= round(2*z2/Tb);
%   if(A_dem1>Ac/2);      % % logic level = (Ac)/2
%     a=1;
%   else(A_dem2>Ac/2);
%     a=0;
%   end
%   y_dem=[y_dem a];
%  
% end
% x_out=y_dem; % output signal;
% 
% % *************** Represent output signal as digital signal ***************
% xx_bit=[];
% for n=1:length(x_out);
%     if x_out(n)==1;
%        xx_bitt=ones(1,nb);
%     else x_out(n)==0;
%         xx_bitt=zeros(1,nb);
%     end
%      xx_bit=[xx_bit xx_bitt];
% end
% t4=Tb/nb:Tb/nb:nb*length(x_out)*(Tb/nb);
% figure(1),subplot(5,2,2)
% plot(t4,xx_bit,'LineWidth',2);grid on;
% axis([ 0 Tb*length(x_out) -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Output signal as digital signal');
% 
% figure(1),subplot(5,2,4)
% plot(t4,xx_bit,'LineWidth',2);grid on;
% axis([ 0 Tb*length(x_out) -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Output signal as digital signal');
% figure(1),subplot(5,2,6)
% plot(t4,xx_bit,'LineWidth',2);grid on;
% axis([ 0 Tb*length(x_out) -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Output signal as digital signal');
% figure(1),subplot(5,2,8)
% plot(t4,xx_bit,'LineWidth',2);grid on;
% axis([ 0 Tb*length(x_out) -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Output signal as digital signal');
% figure(1),subplot(5,2,10)
% plot(t4,xx_bit,'LineWidth',2);grid on;
% axis([ 0 Tb*length(x_out) -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Output signal as digital signal');
% 
% % **************************** end of program *****************************
