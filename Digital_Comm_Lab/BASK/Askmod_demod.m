clc;
clear all;
close all;
% ********************* Define transmitted signal *************************
N=10; %bit sayýsý
% x_inp= round(rand(1,N));  % rastgele 0 ve 1 bitleri 
 x_inp=[ 1 0 0 1 1 0 1 0 0 1]; 

Tb=0.0001; %bit periyodu 
%% input signal
x_bit=[]; 
nb=100; % bbit/bit
for n=1:1:N    
    if x_inp(n)==1;   
       x_bitt=ones(1,nb);
    else x_inp(n)==0;
        x_bitt=zeros(1,nb);
    end
     x_bit=[x_bit x_bitt]; %giris sinyali olusturuldu.
end
t1=Tb/nb:Tb/nb:nb*N*(Tb/nb); % zaman ayarladýk 
f1 = figure(1);
set(f1,'color',[1 1 1]);
subplot(4,1,1);
plot(t1,x_bit,'lineWidth',2);grid on;
axis([ 0 Tb*N -0.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Input signal');
%%  BASK modulasyonu
Ac1=15; % 1 biti için taþýyýcý genliði
Ac2=10; % 0 biti için taþýyýcý genliði
mc=10;  % fc>>fs fc=mc*fs fs=1/Tb 
fc=mc*(1/Tb); % taþýyýcý frekans
t2=Tb/nb:Tb/nb:Tb;                 
t2L=length(t2);
x_mod=[];
for (i=1:1:N)
    if (x_inp(i)==1)
        x_mod0=Ac1*cos(2*pi*fc*t2);%1 için modulasyon iþarti
    else
        x_mod0=Ac2*cos(2*pi*fc*t2);%0 için modulasyon iþareti
    end
    x_mod=[x_mod x_mod0]; %moduleli BASK iþareti
end
t3=Tb/nb:Tb/nb:Tb*N;
% subplot(4,1,2);
% plot(t3,x_mod);
% xlabel('Time(sec)');
% ylabel('Amplitude(volt)');
% title('BASK modulated signal ');
%% AWGN channel model
n1=awgn(x_mod,20);
subplot(4,1,2);
plot(t1,n1,'LineWidth',2);
title('SNR=0 dB')
xlabel('time(sec)');
ylabel('m(t)')

n2=awgn(x_mod,40);
subplot(4,1,3);
plot(t1,n2,'LineWidth',2);
title('SNR=20 dB')
xlabel('time(sec)');
ylabel('m(t)')
%% bu kýsma bak
% NFFT =2^nextpow2(length(x_mod));
%   N=NFFT;
% df=1024000/N;           %normalize frekan deðeri
% ay=-N/2+1:1:N/2;          
% f=ay*df;           % x ekseninde iþaretin frekans deðerlerinin gösterilmesi
%  
% figure
% 
% zey=abs(fftshift(fft(x_mod,NFFT)));
% %stem(f,zey);
% plot(f,zey)


% %********************* iletilen sinyal  x ******************************
% x=x_mod;
% % ********************* Channel model h and w *****************************
% h=1; % zayýflama
% w=0; % gürültü
% % ********************* alýnan sinyal y *********************************
% y=h.*x+w;
% %% BASK demodulasyon
% y_dem=[];
% for n=t2L:t2L:length(y)
%   t=Tb/nb:Tb/nb:Tb;
%   c=cos(2*pi*fc*t); % taþýyýcý frekans 
%   y_dem0=c.*y((n-(t2L-1)):n); %BPF ile tasýyýcý carp
%   t4=Tb/nb:Tb/nb:Tb;
%   z=trapz(t4,y_dem0); % integral ilinti alýcýsý 
%   A_dem=round((2*z/Tb));                                     
%   if(A_dem>((Ac1+Ac2)/2)) % logic level = (Ac1+Ac2)/2
%     A=1;
%   else
%     A=0;
%   end
%   y_dem=[y_dem A];
% end
% x_out=y_dem; % çýkýþ sinyali
% % çýkýþ sinyali
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
% subplot(3,1,3)
% plot(t4,xx_bit,'LineWidth',2);grid on;
% axis([ 0 Tb*length(x_out) -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Demodulated signal');
% % **************************** end of program *****************************
