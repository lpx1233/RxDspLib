clear all;
close all;
clc;

OriginalSignal = importdata('.\Original Data\Original_Data.txt');

OriginalData_port1 = OriginalSignal;
shiftnum = 6;
% OriginalData_port2 = [ones(shiftnum, 1); ~(OriginalSignal(1 : length(OriginalSignal) - shiftnum))];
OriginalData_port2 = [~(OriginalSignal(shiftnum + 1 : end)); ...
											~(OriginalSignal(1 : shiftnum))];
% OriginalData_port2 = 1 - OriginalSignal;
OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
% OriginalData = OriginalData_port1;
SampleRate = 400e9;
OSCRate = 80e9;
DataRate = 25e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = importdata('.\Sampled Data\20170428\woSOA wDSF BtB\-12_6.txt');
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
step(ed1, SampledSignal);

% DownSampledData = SampledSignal(1 : OverSamplingRatio : end, 1);
%
% CorrelationResult = conv(DownSampledData(1:end), conj(OriginalData(end:-1:1)));
% [MaxCorr, index] = max(abs(CorrelationResult));
% plot(CorrelationResult)
% find(abs(CorrelationResult) > 2500)
