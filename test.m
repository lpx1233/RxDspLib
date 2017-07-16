clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

RawSignal = importdata('.\Sampled Data\RoF\wireless\rof_btb_eml20170713.dat');

SampleRate = 600e9;
OSCRate = 120e9;
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = 100 * (RawSignal(:, 2) - mean(RawSignal(:, 2)));
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

t = RawSignal(:, 1);
t = interp1(1:length(t), t, 1:OSCRate/SampleRate:length(t)+1);
t = t';
t = t(1:end-1);

index = find(abs(t) < 1e-9);

fc = 25.5e9;
lo = cos(2*pi*fc*t) + i * -sin(2*pi*fc*t);
% plot(t(index), SampledSignal(index), t(index), real(lo(index)), 'r')

r = lowPassFilter25G(SampledSignal .* lo);

figure;
subplot(2,1,1);
plot(t(index), real(r(index)))
subplot(2,1,2);
plot(t(index), imag(r(index)))

MatchFilter = ones(SampleRate/DataRate, 1);
r = conv(r, MatchFilter, 'same');

figure;
subplot(2,1,1);
plot(t(index), real(r(index)))
subplot(2,1,2);
plot(t(index), imag(r(index)))

scatterplot(r)

carsync = comm.CarrierSynchronizer('Modulation', 'PAM', 'SamplesPerSymbol', OverSamplingRatio);
symsync = comm.SymbolSynchronizer(...
    'TimingErrorDetector', 'Early-Late (non-data-aided)', ...
    'SamplesPerSymbol', OverSamplingRatio);

r = step(carsync, r);

scatterplot(r)

% r = step(symsync, r);
% step(cd3, r)

OriginalSignal = importdata('.\Original Data\Original_Data.txt');
OriginalData_port1 = OriginalSignal;
% shiftnum = 67;
% OriginalData_port2 = [~(OriginalSignal(shiftnum + 1 : end)); ...
% 											~(OriginalSignal(1 : shiftnum))];
shiftnum = 58;
OriginalData_port2 = [~(OriginalSignal(end - shiftnum + 1 : end));
											~(OriginalSignal(1 : end - shiftnum))];
OriginalData = OriginalData_port1 + 2 * OriginalData_port2;
% OriginalData = OriginalData_port1;
SampledSignal = -real(r);
CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
% Need 40min to extract signal
parfor i = 1 : length(CorrelationResult)
  fprintf('Generating %dth result\n', i);
  CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
end
plot(CorrelationResult)
