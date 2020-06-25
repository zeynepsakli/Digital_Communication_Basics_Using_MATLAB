clear all;
close all;
fs=10000; %tabanbannt modüle iþaretin örnekleme frekansý
fd=100; %mod. öncesi ikili iþaretin örnekleme frekansý
b=[0 0 0 0 0 0 1 1 0 0 0 0 0 0]; %binary bit dizini
bnrz=tekkutuplu(b,fd,fs,'nrz'); %nrz iþareti
subplot(421)
plot(bnrz,'b');
stairs(bnrz,'LineWidth',2); %ikili bilgi dizinini çizdir
title('Zamanda kare darbe');
axis([0 1000 0 1.5]);
xlabel('t');
ylabel('amplitude');
fx=fft(bnrz,13000);     %frekans spektrumu hesapla
 subplot(422)                         %130. örnek fb=1/tb oranýna karþýlýk gelir.
plot(abs(fftshift(fx)));   %genlik yanýtý
title('Darbe iþaretinin genlik spektrumu');
xlabel('frekans Hz');
ylabel('amplitude');


fr=fx;    %alýnaný iletilene esitledik
fr(130:13000-130)=0; %kanalýn bantsýnýrlý olmasýný alýnan spektruma yansýt
subplot(423)
plot(abs(fftshift(fr)),'b');
title('Bant sýnýrlý kanalýn sonunda alýnan iþaretin genlik spektrumu')
xlabel('frekans Hz');
ylabel('amplitude');
a=real(ifft(fr,13000)); %alýnan isareti zaman uzayýnda olustur
a=a(1:length(bnrz)); %alýnan iþaretin bilgi kýsmýný kes 
subplot(424)
plot(a,'LineWidth',2);
title('Bant sýnýrlý kanalýn etkisiyle alýcýda alýnan iþaret biçimi ')
axis([0 1000 0 1.5]);
xlabel('t');
ylabel('amplitude')
