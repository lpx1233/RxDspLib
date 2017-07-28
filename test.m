clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

RawSignal = importdata('.\Sampled Data\RoF\wireless\40km\20170727tdc\C3rof 3dBm00000.dat');

SampleRate = 600e9;
OSCRate = 120e9;
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = (RawSignal(:, 1) - mean(RawSignal(:, 1))) / std(RawSignal(:, 1));
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

t = (0:length(SampledSignal)-1)*(1/SampleRate);
t = t';

index = find(t < 1e-9);

fc = 25e9;
lo = cos(2*pi*fc*t) + i * -sin(2*pi*fc*t);

figure;
plot(t(index), SampledSignal(index), t(index), real(lo(index)), 'r')

X = SampledSignal;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1)

r = bandPassFilter12G38G(SampledSignal);

X = r;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1)

r = r .* lo;

X = real(r);
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1)

r = lowPassFilter25G(r);

figure;
subplot(2,1,1);
plot(t(index), real(r(index)))
subplot(2,1,2);
plot(t(index), imag(r(index)))

X = real(r);
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1)

MatchFilter = ones(OverSamplingRatio, 1);
r = conv(r, MatchFilter, 'same');
figure;
subplot(2,1,1);
plot(t(index), real(r(index)))
subplot(2,1,2);
plot(t(index), imag(r(index)))

scatterplot(r)

X = real(r);
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1)

carsync = comm.CarrierSynchronizer('Modulation', 'PAM', ...
  'SamplesPerSymbol', OverSamplingRatio, ...
  'ModulationPhaseOffset', 'Custom', ...
  'CustomPhaseOffset', 0);
[r, phError] = step(carsync, r);

estFreqOffset = diff(phError)*SampleRate/(2*pi);
rmean = cumsum(estFreqOffset)./(1:length(estFreqOffset))';
figure;
plot(rmean)
xlabel('Symbols')
ylabel('Estimated Frequency Offset (Hz)')
grid
scatterplot(r)

% r = -real(r);
% r = (r - mean(r)) / std(r);

% tic
% symsync = comm.SymbolSynchronizer(...
%   'TimingErrorDetector', 'Gardner (non-data-aided)', ...
%   'SamplesPerSymbol', OverSamplingRatio);
% r = step(symsync, r);
% toc
%
% scatterplot(r)

r = -real(r);
r = (r - mean(r)) / std(r);

tic

% OverSamplingRatio = 1;
OriginalSignal = importdata('.\Original Data\Original_Data.txt');
OriginalSignal = (OriginalSignal - 0.5 ) * 2;
OriginalData_port1 = OriginalSignal;
% shiftnum = 67;
% OriginalData_port2 = [~(OriginalSignal(shiftnum + 1 : end)); ...
% 											~(OriginalSignal(1 : shiftnum))];
shiftnum = 58;
OriginalData_port2 = [-(OriginalSignal(end - shiftnum + 1 : end));
											-(OriginalSignal(1 : end - shiftnum))];
% OriginalData = OriginalData_port1 + 2 * OriginalData_port2;
OriginalData = OriginalData_port1;
CorrelationResult = zeros(length(r) - OverSamplingRatio * length(OriginalData) + 1, 1);
% Need 40min to extract signal
parfor i = 1 : length(CorrelationResult)
  % fprintf('Generating %dth result\n', i);
  CorrelationResult(i) = sum(r(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
end
figure;
plot(CorrelationResult)
% TODO: Remove tic/toc
toc

[a, index] = max(CorrelationResult);
ExtractedSignal = r(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);
[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
fprintf('The signal error before equalization\n');
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);
