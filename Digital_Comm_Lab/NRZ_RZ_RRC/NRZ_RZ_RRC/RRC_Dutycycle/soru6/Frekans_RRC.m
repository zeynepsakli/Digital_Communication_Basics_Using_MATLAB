clear all
close all
clc
fs=10;
t=0:1/fs:fs-2;

load('RRC_filtre1.mat');
x=abs(fftshift(fft(filtre,1024)));
figure
plot([-512:511],x);

hold on
load('RRC_filtre2.mat');
y=abs(fftshift(fft(filtre,1024)));

plot([-512:511],y);

hold on
load('RRC_filtre3.mat');
z=abs(fftshift(fft(filtre,1024)));

plot([-512:511],z);
legend('Filtre 1','Filtre 2','Filtre 3');
title('RRC filtreler frekans uzayý');
xlabel('t');
ylabel('amplitude');




