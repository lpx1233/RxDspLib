clear all;
close all;
clc;

OriginalSignal = importdata('.\Original Data\Original_Data.txt');
% OriginalData = 1 - OriginalSignal;
% OriginalData = circshift(OriginalData, 10);
% OriginalData = (OriginalData - 0.5) * 2;

OriginalData_port1 = OriginalSignal;
shiftnum = 6;
% OriginalData_port2 = [ones(shiftnum, 1); ~(OriginalSignal(1 : length(OriginalSignal) - shiftnum))];
OriginalData_port2 = [~(OriginalSignal(end - shiftnum + 1 : end));
											~(OriginalSignal(1 : end - shiftnum))];
% OriginalData_port2 = 1 - OriginalSignal;
OriginalData = OriginalData_port1 + 2 * OriginalData_port2;
% OriginalData = OriginalData_port1;
SampleRate = 400e9;
OSCRate = 80e9;
DataRate = 25e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = importdata('.\Sampled Data\day2\10dml25km10Apd\-16dsf.txt');
% eyediagram(SampledSignal(1:100000), 4*OverSamplingRatio, 2*OverSamplingRatio, 0.5*OverSamplingRatio);
% grid on;

SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
DownSampledData = SampledSignal(1 : OverSamplingRatio : end, 1);

CorrelationResult = conv(DownSampledData(1:end), conj(OriginalData(end:-1:1)));
[MaxCorr, index] = max(abs(CorrelationResult));
plot(CorrelationResult)
% find(abs(CorrelationResult) > 2500)
