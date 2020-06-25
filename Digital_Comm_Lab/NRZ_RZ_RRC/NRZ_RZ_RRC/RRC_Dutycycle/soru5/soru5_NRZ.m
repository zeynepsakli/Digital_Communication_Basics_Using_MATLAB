clear all
close all
clc

Nsym = 6;           % Filter span in symbol durations
beta = 1;         % Roll-off factor
sampsPerSym = 10;   % Upsampling factor

rctFilt = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Normal', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'OutputSamplesPerSymbol', sampsPerSym)

% Normalize to obtain maximum filter tap value of 1
b = coeffs(rctFilt);
rctFilt.Gain = 1/max(b.Numerator);

% Visualize the impulse response
%fvtool(rctFilt, 'Analysis', 'impulse')

% Parameters
DataL = 20;             % Data length in symbols
R = 1000;               % Data rate
Fs = R * sampsPerSym;   % Sampling frequency

% Create a local random stream to be used by random number generators for
% repeatability
hStr = RandStream('mt19937ar', 'Seed', 0);

% Generate random data
x = 2*randi(hStr, [0 1], DataL, 1)-1;
% Time vector sampled at symbol rate in milliseconds
tx = 1000 * (0: DataL - 1) / R;

% Filter
yo = rctFilt([x; zeros(Nsym/2,1)]);
% Time vector sampled at sampling frequency in milliseconds
to = 5000* (0: (DataL+Nsym/2)*sampsPerSym - 1) / Fs;

% Filter group delay, since raised cosine filter is linear phase and
% symmetric.
fltDelay = Nsym / (2*R);
% Correct for propagation delay by removing filter transients
yo = yo(fltDelay*Fs+1:end);
to = 1000 * (0: DataL*sampsPerSym - 1) / Fs;

NRZ_out=[];
Vp=5;

for index=1:size(x,1)
 if x(index)==1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*Vp];
 elseif x(index)==-1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*(-Vp)];
 end
end
subplot(421)
plot(to, NRZ_out, 'r-'); hold off;
title('Kutuplu NRZ transmit data')
% Set axes and labels.
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend('Transmitted Data', 'NRZ Coded Data', 'Unipolar Line Coded Data');
fpass=500; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
subplot(423)
plot(to, NRZ_channel_out, 'b'); hold on;
title('Kutuplu NRZ kanal çýkýþý fc=500 Hz rb=1kbps')
%stem(tx, x, 'kx'); hold off;
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');
fpass=1000; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
subplot(425)
plot(to, NRZ_channel_out, 'b'); hold on;
title('Kutuplu NRZ kanal çýkýþý fc=1000 Hz rb=1kbps')
%stem(tx, x, 'kx'); hold off;
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');
fpass=2000; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
subplot(427)
plot(to, NRZ_channel_out, 'b'); hold on;
title('Kutuplu NRZ kanal çýkýþý fc=2000 Hz rb=1kbps')
%stem(tx, x, 'kx'); hold off;
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');
%% tek kutup icin
Nsym = 6;           % Filter span in symbol durations
beta = 0.1;         % Roll-off factor
sampsPerSym = 10;   % Upsampling factor

rctFilt = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Normal', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'OutputSamplesPerSymbol', sampsPerSym)

% Normalize to obtain maximum filter tap value of 1
b = coeffs(rctFilt);
rctFilt.Gain = 1/max(b.Numerator);

% Visualize the impulse response
%fvtool(rctFilt, 'Analysis', 'impulse')

% Parameters
DataL = 20;             % Data length in symbols
R = 1000;               % Data rate
Fs = R * sampsPerSym;   % Sampling frequency

% Create a local random stream to be used by random number generators for
% repeatability
hStr = RandStream('mt19937ar', 'Seed', 0);

% Generate random data
x = 2*randi(hStr, [0 1], DataL, 1)-1;
% Time vector sampled at symbol rate in milliseconds
tx = 1000 * (0: DataL - 1) / R;

% Filter
yo = rctFilt([x; zeros(Nsym/2,1)]);
% Time vector sampled at sampling frequency in milliseconds
to = 5000* (0: (DataL+Nsym/2)*sampsPerSym - 1) / Fs;

% Filter group delay, since raised cosine filter is linear phase and
% symmetric.
fltDelay = Nsym / (2*R);
% Correct for propagation delay by removing filter transients
yo = yo(fltDelay*Fs+1:end);
to = 1000 * (0: DataL*sampsPerSym - 1) / Fs;

NRZ_out=[];
Vp=5;

%% Tek kutuplu NRZ
for index=1:size(x,1)
 if x(index)==1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*Vp];
 elseif x(index)==-1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*0];
 end
end
subplot(422)
plot(to, NRZ_out, 'r-'); hold off;
title('Tek Kutuplu NRZ transmit data')
% Set axes and labels.
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend('Transmitted Data', 'NRZ Coded Data', 'Unipolar Line Coded Data');
fpass=500; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
subplot(424)
plot(to, NRZ_channel_out, 'b'); hold on;
title(' Tek Kutuplu NRZ kanal çýkýþý fc=500 Hz rb=1kbps')
%stem(tx, x, 'kx'); hold off;
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');
fpass=1000; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
subplot(426)
plot(to, NRZ_channel_out, 'b'); hold on;
title('Tek Kutuplu NRZ kanal çýkýþý fc=1000 Hz rb=1kbps')
% stem(tx, x, 'kx'); hold off;
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');
fpass=2000; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
subplot(428)
plot(to, NRZ_channel_out, 'b'); hold on;
title('Tek Kutuplu NRZ kanal çýkýþý fc=2000 Hz rb=1kbps')
% stem(tx, x, 'kx'); hold off;
axis([0 25 -6 6]);  xlabel('Time (ms)'); ylabel('Amplitude');

