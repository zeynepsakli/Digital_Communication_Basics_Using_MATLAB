clear all
close all
clc

Fs=1500; 
t=0:1/Fs:1-1/Fs;
yt=(cos(2*pi*10*t)+2*cos(2*pi*30*t));

NF=2048;     %%Frekans ayarý için
 df=-NF/2:1:NF/2-1;
y=Fs/NF; 
fY=y*df;

%% 20
impulseTrain=zeros(size(t));
impulseTrain(1:Fs/20:end)=1; %fö=60 dürtü katarý 
Sampled2=yt.*impulseTrain; %ideal örnekleme için çarpýyoruz

frequency2=fftshift(fft(impulseTrain,NF)); %fft islemi
y2_mag=abs(frequency2); %abs deðeri ile çizdirmek icin
subplot(321)
plot(fY,y2_mag);
title('Ýdeal için Örnekleme iþareti 20Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')
%% 60
impulseTrain=zeros(size(t));
impulseTrain(1:Fs/60:end)=1; %fö=60 dürtü katarý 
Sampled2=yt.*impulseTrain; %ideal örnekleme için çarpýyoruz

frequency2=fftshift(fft(impulseTrain,NF)); %fft islemi
y2_mag=abs(frequency2); %abs deðeri ile çizdirmek icin
subplot(322)
plot(fY,y2_mag);
title('Ýdeal için Örnekleme iþareti 60Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')
%% 100
impulseTrain=zeros(size(t));
impulseTrain(1:Fs/100:end)=1; %fö=60 dürtü katarý 

frequency2=fftshift(fft(impulseTrain,NF)); %fft islemi
y2_mag=abs(frequency2); %abs deðeri ile çizdirmek icin
subplot(323)
plot(fY,y2_mag);
title('Ýdeal için Örnekleme iþareti 100Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')
%% 200
impulseTrain=zeros(size(t));
impulseTrain(1:Fs/200:end)=1; %fö=60 dürtü katarý 
Sampled2=yt.*impulseTrain; %ideal örnekleme için çarpýyoruz

frequency2=fftshift(fft(impulseTrain,NF)); %fft islemi
y2_mag=abs(frequency2); %abs deðeri ile çizdirmek icin
subplot(324)
plot(fY,y2_mag);
title(' Ýdeal için Örnekleme iþareti 60Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')


