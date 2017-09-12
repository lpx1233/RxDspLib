clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

OSCRate = 80e9;
DataRate = 25e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;
FileName = '.\Sampled Data\50G PAM4\201708\826data\soa-8-10\-15.txt';

SampledSignal = importdata(FileName);
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);

ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
step(ed, SampledSignal);

% Sequence Extraction
tic
OriginalData = importdata('.\Sampled Data\50G PAM4\201708\826data\data sequence.mat');
OriginalData = OriginalData';
OriginalData(find(OriginalData == 1)) = 3;
OriginalData(find(OriginalData == 0.68)) = 1;
OriginalData(find(OriginalData == 0.4)) = -1;
OriginalData(find(OriginalData == 0)) = -3;
CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
parfor i = 1 : length(CorrelationResult)
  CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
end
figure;
plot(CorrelationResult);
toc

% BER counting
[a, index] = max(CorrelationResult);
ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);
[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
fprintf('The signal error before equalization\n');
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% EQ and BER counting
tic
% ChannelLen = 101:10:501;
% alpha = [0.001; 0.0003; 0.0001];
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

ChanLen1st = 121;
ChanLen2nd = 21;
ChanLen3rd = 0;
alpha = 0.003;
[EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalData, 'lms', 10, ChanLen1st, alpha, ChanLen2nd, [], ChanLen3rd);
figure;
plot(costs);
title('Curve of Convergence');
xlabel('Epoch'); ylabel('Cost');

EqualizedSignalUS = resample(EqualizedSignal, OverSamplingRatio, 1);
ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(EqualizedSignalUS), max(EqualizedSignalUS)]);
step(ed1, EqualizedSignalUS);

[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
fprintf('Equalization Setup: Volterra LMS Channel Length Setup is [%d %d %d], alpha is %f\n', ChanLen1st, ChanLen2nd, ChanLen3rd, alpha);
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);
toc
