clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

% TODO: Remove tic/toc
tic

RawSignal = importdata('.\Sampled Data\RoF\wireless\40km\C3rof 3dBm00000.dat');

SampleRate = 600e9;
OSCRate = 120e9;
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = (RawSignal(:, 2) - mean(RawSignal(:, 2))) / std(RawSignal(:, 2));
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

t = RawSignal(:, 1);
t = interp1(1:length(t), t, 1 : OSCRate/SampleRate : length(t)+1);
t = t';
t = t(1:end-1);

index = find(abs(t) < 1e-9);

fc = 25e9;
lo = cos(2*pi*fc*t) + i * -sin(2*pi*fc*t);
% figure;
% plot(t(index), SampledSignal(index), t(index), real(lo(index)), 'r')

r = lowPassFilter25G(SampledSignal .* lo);
% figure;
% subplot(2,1,1);
% plot(t(index), real(r(index)))
% subplot(2,1,2);
% plot(t(index), imag(r(index)))

MatchFilter = ones(OverSamplingRatio, 1);
r = conv(r, MatchFilter, 'same');
% figure;
% subplot(2,1,1);
% plot(t(index), real(r(index)))
% subplot(2,1,2);
% plot(t(index), imag(r(index)))
% scatterplot(r)

% carsync = comm.CarrierSynchronizer('Modulation', 'PAM', ...
%   'SamplesPerSymbol', OverSamplingRatio, ...
%   'ModulationPhaseOffset', 'Custom', ...
%   'CustomPhaseOffset', 0);
% r = step(carsync, r);
r = r(find(isnan(r) == 0));
% scatterplot(r)

r = real(r);
r = (r - mean(r)) / std(r);
% symsync = comm.SymbolSynchronizer(...
%   'TimingErrorDetector', 'Early-Late (non-data-aided)', ...
%   'SamplesPerSymbol', OverSamplingRatio);
% tic
% r = step(symsync, r);
% toc
% step(cd3, r)

% TODO: Remove tic/toc
toc
tic

OriginalSignal = importdata('.\Original Data\Original_Data.txt');
OriginalSignal = (OriginalSignal - 0.5 ) * 2;
OriginalData_port1 = OriginalSignal;
% shiftnum = 67;
% OriginalData_port2 = [~(OriginalSignal(shiftnum + 1 : end)); ...
% 											~(OriginalSignal(1 : shiftnum))];
shiftnum = 58;
OriginalData_port2 = [~(OriginalSignal(end - shiftnum + 1 : end));
											~(OriginalSignal(1 : end - shiftnum))];
% OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
OriginalData = OriginalData_port1;
CorrelationResult = zeros(length(r) - OverSamplingRatio * length(OriginalData) + 1, 1);
% Need 40min to extract signal
parfor i = 1 : length(CorrelationResult)
  % fprintf('Generating %dth result\n', i);
  CorrelationResult(i) = sum(r(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
end
plot(CorrelationResult)
% prbdet = comm.PreambleDetector('Input', 'Symbol', 'Detections', 'All', ...
%   'Threshold', 100, 'Preamble', OriginalData);
% idx = prbdet(r);
% TODO: Remove tic/toc
toc
