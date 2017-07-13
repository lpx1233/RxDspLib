clear all;
close all;
clc;

OriginalSignal = importdata('.\Original Data\Original_Data.txt');

OriginalData_port1 = OriginalSignal;
% shiftnum = 67;
% OriginalData_port2 = [~(OriginalSignal(shiftnum + 1 : end)); ...
% 											~(OriginalSignal(1 : shiftnum))];

shiftnum = 61;
OriginalData_port2 = [~(OriginalSignal(end - shiftnum + 1 : end));
											~(OriginalSignal(1 : end - shiftnum))];
% OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
OriginalData = OriginalData_port1;

SampleRate = 400e9;
OSCRate = 80e9;
DataRate = 25e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = importdata('.\Sampled Data\50G PAM4\20170428\ebtb.txt');
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
% SampledSignal = SampledSignal - mean(SampledSignal);
% ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
% step(ed1, SampledSignal);

% DownSampledData = SampledSignal(1 : OverSamplingRatio : end, 1);

% CorrelationResult = conv(SampledSignal(1:end), conj(OriginalData(end:-1:1)));
CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
for i = 1 : length(CorrelationResult)
  CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
end
[MaxCorr, index] = max(abs(CorrelationResult));
ExtractedSignalUS = SampledSignal(index-length(OriginalData)*OverSamplingRatio+1 : OverSamplingRatio : index);
plot(CorrelationResult)
% find(abs(CorrelationResult) > 2500)
