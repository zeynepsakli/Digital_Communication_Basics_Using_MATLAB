
clc;
clear all;
close all;
x=[ 1 0 0 1 1 0 1 0 0 1];                                    % Binary Information
bp=0.0001;                                                    % bit period
disp(' Binary information at Trans mitter :');
disp(x);
%XX representation of transmitting binary information as digital signal XXX
bit=[]; 
for n=1:1:length(x)
    if x(n)==1;
       se=ones(1,100);
    else x(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t1=bp/100:bp/100:100*length(x)*(bp/100);
subplot(3,1,1);
plot(t1,bit,'lineWidth',2.5);grid on;
axis([ 0 bp*length(x) -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('Input Signal');
%XXXXXXXXXXXXXXXXXXXXXXX Binary-ASK modulation XXXXXXXXXXXXXXXXXXXXXXXXXXX%
A1=15;                      % Amplitude of carrier signal for information 1
A2=5;                       % Amplitude of carrier signal for information 0
br=1/bp;                                                         % bit rate
f=br*10;                                                 % carrier frequency 
t2=bp/99:bp/99:bp;                 
ss=length(t2);
m=[];
for (i=1:1:length(x))
    if (x(i)==1)
        y=A1*cos(2*pi*f*t2);
    else
        y=A2*cos(2*pi*f*t2);
    end
    m=[m y];
end
t3=bp/99:bp/99:bp*length(x);
fs=1/bp;
NFFT =2^nextpow2(length(m));
N=NFFT;
df=fs/N;           %normalize frekan deðeri
ay=1:N/2;          
fx=(ay*df);           % x ekseninde iþaretin frekans deðerlerinin gösterilmesi



figure
 a_fft=abs(fft(m,1024));
 plot(fx,a_fft(1:N/2));
% subplot(3,1,2);
% plot(t3,m);
% xlabel('time(sec)');
% ylabel('amplitude(volt)');
% title('BASK Modulated Signal');
% zey=fft(m);
% nep=abs(fftshift (zey));
% subplot(313)
% plot(t3,nep);


%XXXXXXXXXXXXXXXXXXXX Binary ASK demodulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% mn=[];
% for n=ss:ss:length(m)
%   t=bp/99:bp/99:bp;
%   y=cos(2*pi*f*t);                                        % carrier siignal 
%   mm=y.*m((n-(ss-1)):n);
%   t4=bp/99:bp/99:bp;
%   z=trapz(t4,mm)                                              % intregation 
%   zz=round((2*z/bp))                                     
%   if(zz>7.5)                                  % logic level = (A1+A2)/2=7.5
%     a=1;
%   else
%     a=0;
%   end
%   mn=[mn a];
% end
% disp(' Binary information at Reciver :');
% disp(mn);
% %XXXXX Representation of binary information as digital signal which achived 
% %after ASK demodulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% bit=[];
% for n=1:length(mn);
%     if mn(n)==1;
%        se=ones(1,100);
%     else mn(n)==0;
%         se=zeros(1,100);
%     end
%      bit=[bit se];
% end
% t4=bp/100:bp/100:100*length(mn)*(bp/100);
% subplot(3,1,2)
% plot(t4,bit,'LineWidth',2.5);grid on;
% axis([ 0 bp*length(mn) -.5 1.5]);
% ylabel('amplitude(volt)');
% xlabel(' time(sec)');
% title('Demodulated Signal');
