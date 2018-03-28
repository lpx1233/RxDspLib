clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

OSCRate = 80e9;
DataRate = 25e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;
FileDir = '.\Sampled Data\201801\O-band DML\40G PD\obtb\';

% TODO import data


% BER counting
[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
fprintf('The signal error before equalization\n');
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% EQ and BER counting
tic
ChannelLen = 151;
alpha = 0.003;
[EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalData, 'lms', ChannelLen, alpha, 5);
% figure;
% plot(costs);
% title('Curve of Convergence');
% xlabel('Epoch'); ylabel('Cost');
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
fprintf('Equalization Setup: Linear LMS Channel Length is %d, alpha is %f\n', ChannelLen, alpha);
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% ChanLen1st = 121;
% ChanLen2nd = 21;
% ChanLen3rd = 0;
% alpha = 0.003;
% [EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalData, 'lms', 10, ChanLen1st, alpha, ChanLen2nd, [], ChanLen3rd);
% figure;
% plot(costs);
% title('Curve of Convergence');
% xlabel('Epoch'); ylabel('Cost');

% EqualizedSignalUS = resample(EqualizedSignal, OverSamplingRatio, 1);
% ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(EqualizedSignalUS), max(EqualizedSignalUS)]);
% step(ed1, EqualizedSignalUS);
%
% [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
% fprintf('Equalization Setup: Volterra LMS Channel Length Setup is [%d %d %d], alpha is %f\n', ChanLen1st, ChanLen2nd, ChanLen3rd, alpha);
% fprintf('Bit error num: %d\n', BitErrorNum);
% fprintf('SER: %e\n', SymErrorRate);
% fprintf('BER: %e\n', BitErrorRate);
toc
