clear all;
close all;
clc;
%Nb is the number of bits to be transmitted
Nb=10000;
% Generate Nb bits randomly
b=rand(1,Nb)>0.5;
%Rb is the bit rate in bits/second
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
% %% kutuplu %50 1000
% for index=1:size(b,2)
%  if b(index)==1
%  RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*Vp]; %50 DBO
%  elseif b(index)==-1
%  RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*(-Vp)];
%  end
% end
% [Pxx f]=pwelch(RZ_out)
% subplot(421)
% plot(f/pi*fs/2,db(abs(Pxx))),grid
% title('PSD-RZ-Kutuplu %50 DBO Rb=1kbps');
% xlabel('Frequency (Hz)')
% %% kutuplu %50 3000
% Rb=3000;
% fs=10*Rb;
% RZ_out=[];
% Vp=5;
% for index=1:size(b,2)
%  if b(index)==1
%  RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*Vp];
%  elseif b(index)==-1
%  RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*(-Vp)];
%  end
% end
% [Pxx f]=pwelch(RZ_out)
% subplot(422)
% plot(f/pi*fs/2,db(abs(Pxx))),grid
% title('RZ-Kutuplu %50 DBO PSD Rb=3kbps');
% xlabel('Frequency (Hz)')
 
%% Tekkutuplu %50 1000
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(421)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('RZ-Tek Kutuplu %50 DBO PSD Rb=1kbps');
xlabel('Frequency (Hz)')
%% Tekkutuplu %50 3000
Rb=3000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(422)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('RZ-Tek Kutuplu %50 DBO PSD Rb=3kbps');
xlabel('Frequency (Hz)')
%% devam
Nb=10000;
% Generate Nb bits randomly
b=rand(1,Nb)>0.5;
%Rb is the bit rate in bits/second
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
% %% kutuplu %20 1000
% for index=1:size(b,2)
%  if b(index)==1
%  RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*Vp];
%  elseif b(index)==-1
%  RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*(-Vp)];
%  end
% end
% [Pxx f]=pwelch(RZ_out)
% subplot(423)
% plot(f/pi*fs/2,db(abs(Pxx))),grid
% title('RZ-Kutuplu %20 DBO PSD Rb=1kbps');
% xlabel('Frequency (Hz)')
% %% kutuplu %20 3000
% Rb=3000;
% fs=10*Rb;
% RZ_out=[];
% Vp=5;
% for index=1:size(b,2)
%  if b(index)==1
%  RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*Vp];
%  elseif b(index)==-1
%  RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*(-Vp)];
%  end
% end
% [Pxx f]=pwelch(RZ_out)
% subplot(424)
% plot(f/pi*fs/2,db(abs(Pxx))),grid
% title('RZ-Kutuplu %20 DBO PSD Rb=3kbps');
% xlabel('Frequency (Hz)')
 
%% Tekkutuplu %20 1000
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(423)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('RZ-Tek Kutuplu %20 DBO PSD Rb=1kbps');
xlabel('Frequency (Hz)')
%% Tekkutuplu %20 3000
Rb=3000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(424)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('RZ-Tek Kutuplu %20 DBO PSD Rb=3kbps');
xlabel('Frequency (Hz)')

