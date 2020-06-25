 close all;
 clear all;
 clc;
 
fs=1500  %Ödev pdf'inde söylendi
tx=0:1/fs:1
y=cos(2*pi*10*tx)+2*cos(2*pi*30*tx);
t = 0 : 1/1500 : 1; % 1 saniye için 1 kHz örneklemne frekansý
impulseTrain=zeros(size(t));
impulseTrain(1:fs/20:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
nx = [0:1/fs:t];
nx = nx(1:end - 1);
fc=20
duty=50
s2 = square(2*pi*20*nx,50); %darbe katarý
s2(find(s2<0)) = 0;
period_samp = length(nx)/fc;
ind = [1:period_samp:length(nx)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(nx));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

L=8; %n=3bit 2^n=8 seviye

t=0:1/1500:1;
xmax=abs(max(pam)); 
adim_araligi=2*xmax/L;
partition=-(L/2-1)*adim_araligi:adim_araligi:(L/2-1)*adim_araligi; %giriþ seviyesi
codebook=-((L-1)*adim_araligi/2):adim_araligi:((L-1)*adim_araligi/2); %çýkýþ seviyesi
[indx xq]=quantiz(pam,partition,codebook); 
%pam sinyalini giriþ ve çýkýþa göre kuantalama islemi
subplot(421)
stem(nx,xq);
title('20 Hz için 3 bit kuantalanmýþ iþaret')
hata=pam-xq;   %hatayý hesaplatan kýsým 
ortalama_hata=sum(hata.^2)./length(hata);   %ortalama kuantalama hatasý
SNR_hesaplat=10*log10(var(pam)/var(hata)) %snr hesabý(Ps/P(n,q))
P=var(pam); %%P(s) görmek icin 
%% 60
impulseTrain(1:fs/60:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
nx = [0:1/fs:t];
nx = nx(1:end - 1);
fc=60
duty=50
s2 = square(2*pi*60*nx,50); %darbe katarý
s2(find(s2<0)) = 0;
period_samp = length(nx)/fc;
ind = [1:period_samp:length(nx)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(nx));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

L=8; %n=3bit 2^n=8 seviye

t=0:1/1500:1;
xmax=abs(max(pam)); %xmax=4;

adim_araligi=2*xmax/L;
partition=-(L/2-1)*adim_araligi:adim_araligi:(L/2-1)*adim_araligi; %giriþ seviyesi

codebook=-((L-1)*adim_araligi/2):adim_araligi:((L-1)*adim_araligi/2); %çýkýþ seviyesi

[indx xq]=quantiz(pam,partition,codebook); 
%pam sinyalini giriþ ve çýkýþa göre kuantalama islemi
subplot(422), stem(nx,xq);
title('60 Hz için 3 bit kuantalanmýþ iþaret')
hata=pam-xq;   %hatayý hesaplatan kýsým
 
ortalama_hata=sum(hata.^2)./length(hata);   %ortalama kuantalama hatasý

SNR_hesaplat=10*log10(var(pam)/var(hata)) %snr hesabý(Ps/P(n,q))

P=var(pam); %%P(s) görmek icin 
%% 100
impulseTrain(1:fs/100:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
nx = [0:1/fs:t];
nx = nx(1:end - 1);
fc=100
duty=50
s2 = square(2*pi*100*nx,50); %darbe katarý
s2(find(s2<0)) = 0;
period_samp = length(nx)/fc;
ind = [1:period_samp:length(nx)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(nx));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

L=8; %n=3bit 2^n=8 seviye
t=0:1/1500:1;
xmax=abs(max(pam)); %xmax=4;

adim_araligi=2*xmax/L;
partition=-(L/2-1)*adim_araligi:adim_araligi:(L/2-1)*adim_araligi; %giriþ seviyesi

codebook=-((L-1)*adim_araligi/2):adim_araligi:((L-1)*adim_araligi/2); %çýkýþ seviyesi

[indx xq]=quantiz(pam,partition,codebook); 
%pam sinyalini giriþ ve çýkýþa göre kuantalama islemi
subplot(423), stem(nx,xq);
title('100 Hz için 3 bit kuantalanmýþ iþaret')
hata=pam-xq;   %hatayý hesaplatan kýsým
 
ortalama_hata=sum(hata.^2)./length(hata);   %ortalama kuantalama hatasý

SNR_hesaplat=10*log10(var(pam)/var(hata)) %snr hesabý(Ps/P(n,q))

P=var(pam); %%P(s) görmek icin 
%% 200
impulseTrain(1:fs/250:end)=1;
Sampled=y.*impulseTrain; %ideal örnekleme için çarpýyoruz
t = 1;  
nx = [0:1/fs:t];
nx = nx(1:end - 1);
fc=250
duty=50
s2 = square(2*pi*250*nx,50); %darbe katarý
s2(find(s2<0)) = 0;
period_samp = length(nx)/fc;
ind = [1:period_samp:length(nx)];
on_samp = ceil(period_samp * 50/100);
pam = zeros(1,length(nx));
for i = 1 : length(ind)
   pam(ind(i):ind(i) + on_samp) = Sampled(ind(i));
end

L=8; %n=3bit 2^n=8 seviye
t=0:1/1500:1;
xmax=abs(max(pam));

adim_araligi=2*xmax/L; %% delta hesaplattýk
partition=-(L/2-1)*adim_araligi:adim_araligi:(L/2-1)*adim_araligi; %giriþ seviyesi

codebook=-((L-1)*adim_araligi/2):adim_araligi:((L-1)*adim_araligi/2); %çýkýþ seviyesi
[indx xq]=quantiz(pam,partition,codebook); 
%pam sinyalini giriþ ve çýkýþa göre kuantalama islemi
subplot(424), stem(nx,xq);
title('200 Hz için 3 bit kuantalanmýþ iþaret')
hata=pam-xq;   %hatayý hesaplatan kýsým
ortalama_hata=sum(hata.^2)./length(hata);   %ortalama kuantalama hatasý
SNR_hesaplat=10*log10(var(pam)/var(hata)) %snr hesabý(Ps/P(n,q))

P=var(pam); %%P(s) görmek icin 


 