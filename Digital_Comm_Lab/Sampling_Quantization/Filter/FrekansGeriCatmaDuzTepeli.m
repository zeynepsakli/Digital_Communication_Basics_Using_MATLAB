clear all
close all
clc
clear all
close all
clc
Fs=1500  %Ödev pdf'inde söylendi
tx=0:1/Fs:1-1/1500
NF=2048;
df=-NF/2:1:NF/2-1;
y=Fs/NF; 
fY=y*df;
y=cos(2*pi*10*tx)+2*cos(2*pi*30*tx) %Orjinal Sinyal
t = 0 : 1/1500 : 1-1/1500; % 1 saniye için 1 kHz örneklemne frekansý
impulseTrain=zeros(size(t));
impulseTrain(1:Fs/20:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
n = [0:1/Fs:1];
n = n(1:end - 1);
fc=20
duty=10
s = square(2*pi*20*n,10);
s(find(s<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * duty/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);

frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(421)
plot(df,y1_mag);
title('Geri Catma Düz Tepeli Örnek %10-DBO fö=20 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')

s2 = square(2*pi*20*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);
frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(422)
plot(df,y1_mag);
title('Geri Catma Düz Tepeli Örnek %50-DBO fö=20 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')

%% 60 Hz

impulseTrain=zeros(size(t));
impulseTrain(1:Fs/60:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
n = [0:1/Fs:1];
n = n(1:end - 1);
fc=60
duty=10
s = square(2*pi*60*n,10);
s(find(s<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * duty/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);

frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(423)
plot(df,y1_mag);
title('Geri Catma Düz Tepeli Örnek %10-DBO fö=60 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')

s2 = square(2*pi*60*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);
frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(424)
plot(df,y1_mag);
title('Geri Catma Düz Tepeli Örnek %50-DBO fö=60 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')

%% 100 Hz

impulseTrain=zeros(size(t));
impulseTrain(1:Fs/100:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
n = [0:1/Fs:1];
n = n(1:end - 1);
fc=100;
duty=10;
s = square(2*pi*100*n,10);
s(find(s<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * duty/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);

frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(425)
plot(df,y1_mag);
title('Geri Catma Düz Tepeli Örnek %10-DBO fö=100 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')

s2 = square(2*pi*100*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);
frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(426)
plot(df,y1_mag);
title('Geri Catma Düz Tepeli Örnek %50-DBO fö=100 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')

%% 200 Hz

impulseTrain=zeros(size(t));
impulseTrain(1:Fs/250:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
n = [0:1/Fs:1];
n = n(1:end - 1);
fc=250;
duty=10;
s = square(2*pi*250*n,10);
s(find(s<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * duty/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);

frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(427)
plot(df,y1_mag);
title('Geri Catma Düz Tepeli Örnek %10-DBO fö=200 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')

s2 = square(2*pi*250*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end
[B,A] = butter(10,30/(1500/2));
%freqz(B,A)
rend=filter(B,A,pam);
frequency1=fftshift(fft(rend,NF));
y1_mag=abs(frequency1);

subplot(428)
plot(df,y1_mag);
title(' Geri Catma Düz Tepeli Örnek %50-DBO fö=200 Hz Frekans Uzayý')
xlabel('frekans Hz')
ylabel('amplitude')


