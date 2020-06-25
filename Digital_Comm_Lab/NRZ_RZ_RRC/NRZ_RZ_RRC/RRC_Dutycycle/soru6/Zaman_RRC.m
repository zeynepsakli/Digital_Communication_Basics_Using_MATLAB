clear all
close all
clc
fs=10; %uysun diye verdim
t=0:1/fs:fs-2; %zaman ayarladým filtre zamanýyla uyuþsun diye
load('RRC_filtre1.mat');
figure
plot(t,filtre);
hold on
load('RRC_filtre2.mat');
plot(t,filtre);
hold on
load('RRC_filtre3.mat');
plot(t,filtre);
legend('filtre 1 azalma fak. 0.5 ','Filtre 2 azalma fak. 1','Filtre 3 azalma fak. 0.25');
title('RRC filtreler zaman uzayý');
xlabel('t');
ylabel('amplitude');





