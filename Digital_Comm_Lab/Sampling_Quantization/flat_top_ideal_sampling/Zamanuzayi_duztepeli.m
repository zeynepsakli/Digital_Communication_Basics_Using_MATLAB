clear all
close all
clc
%%20
fs=1500  %Ödev pdf'inde söylendi
tx=0:1/fs:1
y=cos(2*pi*10*tx)+2*cos(2*pi*30*tx) %Orjinal Sinyal

t = 0 : 1/1500 : 1; % 1 saniye için 1 kHz örneklemne frekansý
impulseTrain=zeros(size(t));
impulseTrain(1:fs/20:end)=1;

Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz

t = 1;  
n = [0:1/fs:t];
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

subplot(421);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %10-DBO fö=20 Hz')
xlabel 't', ylabel Amplitude
s2 = square(2*pi*20*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

subplot(422);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %50-DBO fö=20 Hz')
xlabel 't', ylabel Amplitude
%% 60
impulseTrain=zeros(size(t));
impulseTrain(1:fs/60:end)=1;

Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz

t = 1;  
n = [0:1/fs:t];
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

subplot(423);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %10-DBO fö=60 Hz')
xlabel 't', ylabel Amplitude
s2 = square(2*pi*60*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

subplot(424);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %50-DBO fö=60 Hz')
xlabel 't', ylabel Amplitude

impulseTrain=zeros(size(t));
impulseTrain(1:fs/100:end)=1;

Sampled=y.*impulseTrain;

fc=100
duty=10
s = square(2*pi*100*n,10);
s(find(s<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * duty/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

subplot(425);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %10-DBO fö=100 Hz')
xlabel 't', ylabel Amplitude
s2 = square(2*pi*100*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

subplot(426);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %50-DBO fö=100 Hz')
xlabel 't', ylabel Amplitude

impulseTrain=zeros(size(t));
impulseTrain(1:fs/250:end)=1;

Sampled=y.*impulseTrain;
fc=250
duty=10
s = square(2*pi*250*n,10);
s(find(s<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * duty/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

subplot(427);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %10-DBO fö=200 Hz')
xlabel 't', ylabel Amplitude
s2 = square(2*pi*250*n,50);
s2(find(s2<0)) = 0;
period_samp = length(n)/fc;
ind = [1:period_samp:length(n)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(n));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

subplot(428);plot(n,pam);ylim([-4 4]);
title('Düz Tepeli Örnekleme %50-DBO fö=200 Hz')
xlabel 't', ylabel Amplitude



