OriginalSignal = importdata('.\Original Data\Original_Data.txt');
% OriginalData = 1 - OriginalSignal;
% OriginalData = circshift(OriginalData, 10);
% OriginalData = (OriginalData - 0.5) * 2;

OriginalData_port1 = OriginalSignal;
shiftnum = 7;
% OriginalData_port2 = [ones(shiftnum, 1); ~(OriginalSignal(1 : length(OriginalSignal) - shiftnum))];
OriginalData_port2 = [~(OriginalSignal(shiftnum + 1 : length(OriginalSignal))); ones(shiftnum, 1)];
% OriginalData_port2 = 1 - OriginalSignal;
OriginalData = OriginalData_port1 + 2 * OriginalData_port2;
% OriginalData = OriginalData_port2;
SampleRate = 400e9;
OSCRate = 80e9;
DataRate = 25e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = importdata('.\Sampled Data\50G PAM4 BtB\obtb w dsf -15.txt');

SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
DownSampledData = SampledSignal(1 : OverSamplingRatio : end, 1);

CorrelationResult = conv(DownSampledData(1:end), conj(OriginalData(end:-1:1)));
[MaxCorr, index] = max(abs(CorrelationResult));
plot(CorrelationResult)
% find(abs(CorrelationResult) > 2500)
