clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

RawSignal = importdata('.\Sampled Data\RoF\wireless\40km\20170801_2_4Vpp\C3rof 6-5dBm00000.dat');

SampleRate = 600e9;
OSCRate = 120e9;
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = (RawSignal(:, 1) - mean(RawSignal(:, 1))) / std(RawSignal(:, 1));
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

t = (0:length(SampledSignal)-1)*(1/SampleRate);
t = t';

index = find( (t < 22e-9) & (t > 2e-9));

X = SampledSignal;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,20 * log10(P1))

r = bandPassFilter12G38G(SampledSignal);

figure;
plot(t(index), r(index))

X = r;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,20 * log10(P1))

r = r .^ 2;

X = r;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,20 * log10(P1))

r = lowPassFilter12_5G(r);

figure;
plot(t(index), r(index))

X = r;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = SampleRate;
f = Fs*(0:(L/2))/L;
figure;
plot(f,20 * log10(P1))

r = (r - mean(r)) / std(r);

ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', 48, 'YLimits', [min(r(100000: 200000)), max(r(100000: 200000))]);
step(ed, r(100000: 200000));

tic
% OverSamplingRatio = 1;
OriginalSignal = importdata('.\Original Data\Original_Data.txt');
OriginalSignal = (OriginalSignal - 0.5 ) * 2;
OriginalData_port1 = OriginalSignal;
shiftnum = 58;
OriginalData_port2 = [-(OriginalSignal(end - shiftnum + 1 : end));
											-(OriginalSignal(1 : end - shiftnum))];
OriginalData = OriginalData_port1 + 2 * OriginalData_port2;
% OriginalData = OriginalData_port1;
CorrelationResult = zeros(length(r) - OverSamplingRatio * length(OriginalData) + 1, 1);
parfor i = 1 : length(CorrelationResult)
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

% linear FFE

% ChannelLen = 101:10:501;
% alpha = [0.01; 0.003; 0.001];
% BER = zeros(length(ChannelLen), length(alpha));
% BitError = zeros(length(ChannelLen), length(alpha));
% for i = 1 : length(ChannelLen)
%   for j = 1 : length(alpha)
%     [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalData, 'lms', ChannelLen(i), alpha(j), 5);
%     % figure;
%     % plot(costs);
%     % title('Curve of Convergence');
%     % xlabel('Epoch'); ylabel('Cost');
%
%     [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
%     fprintf('Equalization Setup: Linear LMS Channel Length is %d, alpha is %f\n', ChannelLen(i), alpha(j));
%     fprintf('Bit error num: %d\n', BitErrorNum);
%     fprintf('SER: %e\n', SymErrorRate);
%     fprintf('BER: %e\n', BitErrorRate);
%     BitError(i, j) = BitErrorNum;
%     BER(i, j) = BitErrorRate;
%   end
% end
% meshz(log10(alpha), ChannelLen, -log10(BER));
% title('BER vs ChannelLen & alpha')
% xlabel('log10(alpha)') % x-axis label
% ylabel('ChannelLen') % y-axis label
% zlabel('-log10(BER)') % y-axis label

% tic
%
% ChanLen1st = 301;
% ChanLen2nd = 33;
% ChanLen3rd = 11;
% alpha = 0.003;
%
% [EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalData, 'lms', 10, ChanLen1st, alpha, ChanLen2nd, [], ChanLen3rd);
% figure;
% plot(costs);
% title('Curve of Convergence');
% xlabel('Epoch'); ylabel('Cost');
%
% EqualizedSignalUS = resample(EqualizedSignal, 48, 1);
% ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', 48, 'YLimits', [min(EqualizedSignalUS), max(EqualizedSignalUS)]);
% step(ed1, EqualizedSignalUS);
%
% [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
% fprintf('Equalization Setup: Volterra LMS Channel Length Setup is [%d %d %d], alpha is %f\n', ChanLen1st, ChanLen2nd, ChanLen3rd, alpha);
% fprintf('Bit error num: %d\n', BitErrorNum);
% fprintf('SER: %e\n', SymErrorRate);
% fprintf('BER: %e\n', BitErrorRate);
% toc
