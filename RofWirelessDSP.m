clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

%% Generate original data
OriginalData = generateData();

% Define parameters
SampleRate = 400e9;
OSCRate = 80e9;
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

% Import data
PathToSampledData = '.\Sampled Data\rof_data_Btb.csv';
x = importdata(PathToSampledData);
SampledSignal = x(:, 2);
%% Resample to 400GSps
Time = -10e-6 : (20e-6 / (8e6 - 1)) : 10e-6;
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

% Mix signal with 25GHz lo

% Signal Synchronization and Extraction
[ExtractedSignal, OriginalSignal] = syncAndExtractSignal(SampledSignal, OriginalData, OverSamplingRatio);

[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] = decisionAndCalcBerPAM4(ExtractedSignal, OriginalSignal);
fprintf('The signal error before equalization for %s\n', PathToSampledData);
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% Equalization & BER Counting
% [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 21, 0.01, 5);
% % figure;
% % plot(costs);
% % title('Curve of Convergence');
% % xlabel('Epoch'); ylabel('Cost');
% [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
% fprintf('Equalization Setup: Linear LMS Channel Length is %d, alpha is 0.01\n', 21);
% fprintf('Bit number num: %d \n', BitErrorNum);
% fprintf('SER: %e\n', SymErrorRate);
% fprintf('BER: %e\n', BitErrorRate);
