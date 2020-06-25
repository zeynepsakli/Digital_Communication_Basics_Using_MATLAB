clc;
clear all;
close all;
% ********************* Define transmitted signal *************************
N=10; %bit sayýsý
%x_inp= round(rand(1,N));  % rastgele 0 ve 1 bitleri 
x_inp=[1 0 0 1 1 0 1 0 0 1]; 

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
% f1 = figure(1);
% set(f1,'color',[1 1 1]);
% subplot(4,1,1);
% plot(t1,x_bit,'lineWidth',2);grid on;
% axis([ 0 Tb*N -0.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Input signal');
%%  BASK modulasyonu
Ac1=15; % 1 biti için taþýyýcý genliði
Ac2=10; % 0 biti için taþýyýcý genliði
mc=10;  % fc>>fs fc=mc*fs fs=1/Tb 
fc=mc*(1/Tb); % taþýyýcý frekans
t2=Tb/nb:Tb/nb:Tb;                 
t2L=length(t2);
 x_mod=[];
y_mod=[];
rend_mod=[];
y_dem=[];


for (i=1:1:N)
    if (x_inp(i)==1)
        x_mod0=Ac1*cos(2*pi*fc*t2);%1 için modulasyon iþarti
        y = abs(x_mod0); %Rectifier
        [B,A] = butter(10,fc/1000000,'low'); %LPF için
        rend=filter(B,A,y); %LPF filtre
         else
        x_mod0=Ac2*cos(2*pi*fc*t2);%0 için modulasyon iþareti
        y = abs(x_mod0);
        [B,A] = butter(10,fc/1000000,'low');
        rend=filter(B,A,y);
    end
    %*************COMPARATOR*****************
        if(rend<=9) % threshold verdik zarf dedektorune gore
     A=0;
   else
     A=1;
        end
   
   y_dem=[y_dem A];
   x_out=y_dem; % çýkýþ sinyali

xx_bit=[];
 for n=1:length(x_out);
     if x_out(n)==1;
        xx_bitt=ones(1,nb);
     else x_out(n)==0;
        xx_bitt=zeros(1,nb);
    end
     xx_bit=[xx_bit xx_bitt]; %Comparator  çýkýþý

 end 
    x_mod=[x_mod x_mod0]; %moduleli BASK iþareti
    y_mod=[y_mod y]; %rectifier için
    rend_mod=[rend_mod rend]; %LPF filtre 
end
t3=Tb/nb:Tb/nb:Tb*N;
subplot(4,1,1);
plot(t3,x_mod);
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
title('BASK modulated signal ');

subplot(412)
 plot(t3,y_mod);
grid
title('doðrultucu çýkýþý')
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
subplot(413)
 plot(t3,rend_mod,'LineWidth',1,'Color','r');
 title('LPF çýkýþý')
 axis([0 0.001 0 15]);
 xlabel('Time(sec)');
ylabel('Amplitude(volt)');
grid

t4=Tb/nb:Tb/nb:nb*length(x_out)*(Tb/nb);
subplot(414)
plot(t4,xx_bit,'LineWidth',2,'Color','r');grid on;
axis([ 0 Tb*length(x_out) -0.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Comparator Çýkýþý');



