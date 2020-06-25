clear all
close all
clc
%% x(t)
fa=1500  %Ödev pdf'inde söylendi
tx=0:1/fa:1 %periuyod ayarladým
y=cos(2*pi*10*tx)+2*cos(2*pi*30*tx) %Orjinal Sinyal
subplot(321)
plot(tx,y)
xlabel('t')   %x ekseni
ylabel('amplitude') %y ekseni
title('Sinyal') %baslik 

%% fö=20Hz 
% dürtü katarý oluþturuyoruz.
t = 0 : 1/1500 : 1; % 1 saniye için 1 kHz örneklemne frekansý
impulseTrain=zeros(size(t)); %ÖDevde söylendi 
impulseTrain(1:fa/20:end)=1; %fö=20 olarak dürtü 
xlabel 'Time (s)', ylabel Waveform
subplot(322)
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
plot(t,Sampled)
hold on
plot(t,Sampled,'*') %noktalar koysun diye
title('Ýdeal Örnekleme fö=20 Hz')
%% fö=60 Hz
%t = 0 : 1/1500 : 1; % 1 saniye için 1 kHz örneklemne frekansý
impulseTrain=zeros(size(t));
impulseTrain(1:fa/60:end)=1;
xlabel 'Time (s)', ylabel Waveform
%%pulse train oluþtu. Þimdi kendi sinyalimizi oluþturuyoruz.
y=cos(2*pi*10*t)+2*cos(2*pi*30*t); % 5 Hz frekanslý. Sn’de 5 periyot oluþturacak.
subplot(323)
Sampled2=y.*impulseTrain;
plot(t,Sampled2)
hold on
plot(t,Sampled2,'*')
title('Ýdeal Örnekleme fö=60 Hz')
%% fö=100 Hz
impulseTrain=zeros(size(t));
impulseTrain(1:fa/100:end)=1;
xlabel 'Time (s)', ylabel Waveform
%%pulse train oluþtu. Þimdi kendi sinyalimizi oluþturuyoruz.
y=cos(2*pi*10*t)+2*cos(2*pi*30*t); % 5 Hz frekanslý. Sn’de 5 periyot oluþturacak.
subplot(324)
Sampled2=y.*impulseTrain;
plot(t,Sampled2)
hold on
plot(t,Sampled2,'*')
title('Ýdeal Örnekleme fö=100 Hz')
%% fö=200 Hz
impulseTrain=zeros(size(t));
impulseTrain(1:fa/200:end)=1;
xlabel 'Time (s)', ylabel Waveform
%%pulse train oluþtu. Þimdi kendi sinyalimizi oluþturuyoruz.
y=cos(2*pi*10*t)+2*cos(2*pi*30*t); % 5 Hz frekanslý. Sn’de 5 periyot oluþturacak.
subplot(325)
Sampled2=y.*impulseTrain;
plot(t,Sampled2)
hold on
plot(t,Sampled2,'*')
title('Ýdeal Örnekleme fö=200 Hz')
xlabel 'Time (s)', ylabel Waveform



