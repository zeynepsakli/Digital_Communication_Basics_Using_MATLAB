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
to = 1000* (0: (DataL+Nsym/2)*sampsPerSym - 1) / Fs;

% Filter group delay, since raised cosine filter is linear phase and
% symmetric.
fltDelay = Nsym / (2*R);
% Correct for propagation delay by removing filter transients
yo = yo(fltDelay*Fs+1:end);
to = 1000 * (0: DataL*sampsPerSym - 1) / Fs;
% Plot data.
%stem(tx, x, 'kx'); hold on;
% Plot filtered data.
plot(to, yo, 'r'); hold on;


%Line Coding start
NRZ_out=[];
Vp=1;

for index=1:size(x,1)
 if x(index)==1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*Vp];
 elseif x(index)==-1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*(-Vp)];
 end
end
%Line Coding stop


plot(to, NRZ_out, 'b'); hold off;
% Set axes and labels.
axis([0 25 -2.5 2.5]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend( 'NRZ coded data', 'Unipolar Line Coded Data','Location', 'southeast')



% CHANNEL MODEL
fpass=1000; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
% CHANNEL MODEL

% CHANNEL MODEL
NRZ_channel_out_rc = lowpass(yo,fpass,Fs);
% CHANNEL MODEL


figure
plot(to, NRZ_channel_out, 'r'); hold on;
plot(to, NRZ_channel_out_rc, 'b'); hold on;
%stem(tx, x, 'kx'); hold off;
% Set axes and labels.
axis([0 25 -2.5 2.5]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend('nrz channel out', 'RC coded data', 'Location', 'southeast')
title('Channel Output Signal fc=1000 Hz')
%% 500
fpass=500; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
% CHANNEL MODEL

% CHANNEL MODEL
NRZ_channel_out_rc = lowpass(yo,fpass,Fs);
% CHANNEL MODEL


figure
plot(to, NRZ_channel_out, 'r'); hold on;
plot(to, NRZ_channel_out_rc, 'b'); hold on;
%stem(tx, x, 'kx'); hold off;
% Set axes and labels.
axis([0 25 -2.5 2.5]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend('nrz channel out', 'RC coded data', 'Location', 'southeast')
title('Channel Output Signal fc=500 Hz')

%% 100
fpass=100; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
% CHANNEL MODEL

% CHANNEL MODEL
NRZ_channel_out_rc = lowpass(yo,fpass,Fs);
% CHANNEL MODEL


figure
plot(to, NRZ_channel_out, 'r'); hold on;
plot(to, NRZ_channel_out_rc, 'b'); hold on;
%stem(tx, x, 'kx'); hold off;
% Set axes and labels.
axis([0 25 -2.5 2.5]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend('nrz channel out', 'RC coded data', 'Location', 'southeast')
title('Channel Output Signal fc=100 Hz')
%% 10
fpass=10; %change!
NRZ_channel_out = lowpass(NRZ_out,fpass,Fs);
% CHANNEL MODEL

% CHANNEL MODEL
NRZ_channel_out_rc = lowpass(yo,fpass,Fs);
% CHANNEL MODEL


figure
plot(to, NRZ_channel_out, 'r'); hold on;
plot(to, NRZ_channel_out_rc, 'b'); hold on;
%stem(tx, x, 'kx'); hold off;
% Set axes and labels.
axis([0 25 -2.5 2.5]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend('nrz channel out', 'RC coded data', 'Location', 'southeast')
title('Channel Output Signal fc=10 Hz')




