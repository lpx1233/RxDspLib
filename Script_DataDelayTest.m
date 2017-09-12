clear all;
close all;
clc;

OSCRate = 80e9;
DataRate = 28e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = importdata('.\Sampled Data\50G PAM4\201708\cal\ebtb\C2ebtb00000.dat');
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);

ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
step(ed, SampledSignal);

tic
OriginalSignal = importdata('.\Original Data\prbs15pam4.mat');
OriginalSignal = OriginalSignal';
OriginalSignal = (OriginalSignal - mean(OriginalSignal)) * 6;

OriginalData_port1 = OriginalSignal;
shiftnum = 61;
% OriginalData_port2 = [~(OriginalSignal(shiftnum + 1 : end)); ...
% 											~(OriginalSignal(1 : shiftnum))];
OriginalData_port2 = [~(OriginalSignal(end - shiftnum + 1 : end));
											~(OriginalSignal(1 : end - shiftnum))];
% OriginalData = 2 * OriginalData_port1 + OriginalData_port2;

OriginalData = OriginalSignal;
CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
parfor i = 1 : length(CorrelationResult)
  CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
end
plot(CorrelationResult)
toc
