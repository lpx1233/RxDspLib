clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

%% Import sampled data from DSO
% Defining parameters
SampleRate = 200e9;
OSCRate = 40e9;
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = importdata('.\Sampled Data\RoF\wired\20170802\C1-1800000.dat');
SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
step(ed, SampledSignal);

tic
% OverSamplingRatio = 1;
OriginalSignal = importdata('.\Original Data\Original_Data.txt');
OriginalSignal = (OriginalSignal - 0.5 ) * 2;
OriginalData_port1 = OriginalSignal;
shiftnum = 58;
% OriginalData_port2 = [-(OriginalSignal(end - shiftnum + 1 : end));
%                       -(OriginalSignal(1 : end - shiftnum))];
OriginalData_port2 = [-(OriginalSignal(shiftnum + 1 : end));
                      -(OriginalSignal(1 : shiftnum))];
OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
% OriginalData = OriginalData_port1;
CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
parfor i = 1 : length(CorrelationResult)
  CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
end
figure;
plot(CorrelationResult)
% TODO: Remove tic/toc
toc

[a, index] = max(CorrelationResult);
ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);
[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
fprintf('The signal error before equalization\n');
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

[EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalData, 'lms', 21, 1e-3, 10);
figure;
plot(costs);
title('Curve of Convergence');
xlabel('Epoch'); ylabel('Cost');

EqualizedSignalUS = resample(EqualizedSignal, OverSamplingRatio, 1);
ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(EqualizedSignalUS), max(EqualizedSignalUS)]);
step(ed1, EqualizedSignalUS);

[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
fprintf('Equalization Setup: Linear LMS Channel Length is 21, alpha is 1e-3\n');
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% linear FFE
% ChannelLen = 11:10:201;
% alpha = [0.003; 0.001; 0.0003];
% BER = zeros(length(ChannelLen), length(alpha));
% BitError = zeros(length(ChannelLen), length(alpha));
% for i = 1 : length(ChannelLen)
%   for j = 1 : length(alpha)
%     [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalData, 'lms', ChannelLen(i), alpha(j), 10);
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


% PathToSampledDataFolder = '.\Sampled Data\RoF\wired\20170802\';
% FileList = dir(PathToSampledDataFolder);
% BERwEQ = zeros(length(FileList) - 2, 1);
% BERwoEQ = BERwEQ;
% for i = 1 : length(FileList) - 2
%   PathToSampledData = strcat(PathToSampledDataFolder, FileList(i+2).name);
%   SampledSignal = importdata(PathToSampledData);
%   SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);
%   SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
%   tic
%   OriginalSignal = importdata('.\Original Data\Original_Data.txt');
%   OriginalSignal = (OriginalSignal - 0.5 ) * 2;
%   OriginalData_port1 = OriginalSignal;
%   shiftnum = 58;
%   OriginalData_port2 = [-(OriginalSignal(shiftnum + 1 : end));
%                         -(OriginalSignal(1 : shiftnum))];
%   OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
%   % OriginalData = OriginalData_port1;
%   CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
%   parfor i = 1 : length(CorrelationResult)
%     CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
%   end
%   figure;
%   plot(CorrelationResult)
%   % TODO: Remove tic/toc
%   toc
%
%   [a, index] = max(CorrelationResult);
%   ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);
%   [BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
%   fprintf('The signal error before equalization\n');
%   fprintf('Bit error num: %d\n', BitErrorNum);
%   fprintf('SER: %e\n', SymErrorRate);
%   fprintf('BER: %e\n', BitErrorRate);
%   BERwoEQ(i) = BitErrorRate;
%
% 	[EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalData, 'lms', 21, 1e-3, 10);
% 	% figure;
% 	% plot(costs);
% 	% title('Curve of Convergence');
% 	% xlabel('Epoch'); ylabel('Cost');
%
% 	[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
% 	fprintf('Equalization Setup: Linear LMS Channel Length is 21, alpha is 1e-3\n');
% 	fprintf('Bit number num: %d \n', BitErrorNum);
% 	fprintf('SER: %e\n', SymErrorRate);
% 	fprintf('BER: %e\n', BitErrorRate);
% 	BERwEQ(i) = BitErrorRate;
% end
%
% x = [-12; -13; -14; -15; -17; -18; -19; -20; -21; -22];
% plot(x, log10(BERwEQ), x, log10(BERwoEQ), 'r')
% title('BER curve w(b)/wo(r) FFE');
% xlabel('ROP/dBm'); ylabel('log10(BER)');
