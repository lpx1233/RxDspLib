%% This example shows the work flow for a receiver side dsp: data generation,
%% synchronization, signal extraction, PAM4 decision and BER calculation.
%% Also the eyediagram is drawn, but since the equalizer is 1 sample/sym,
%% the equalized signal eyediagram cannot be drawn. What's more, the error
%% performance of signal before and after equalization is compared.

clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

%% Generate original data
OriginalData = generateData();

%% Import sampled data from DSO
% Defining parameters
SampleRate = 400e9;
OSCRate = 80e9;
DataRate = 25e9;
OverSamplingRatio = SampleRate / DataRate;
% importing and eyediagram drawing
SampledSignal = importdata('.\Sampled Data\20170428\ebtb.txt');
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
eyediagram(SampledSignal(1:100000), 4*OverSamplingRatio, 2*OverSamplingRatio, 0.5*OverSamplingRatio);
grid on;

%% Signal Synchronization and Extraction
[ExtractedSignal, OriginalSignal] = syncAndExtractSignal(SampledSignal, OriginalData, OverSamplingRatio);

%% LMS Equalization
% 101-tap FFE and training for 5 epochs
% [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 1001, 0.001, 5);
% [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalSignal, 'rls', 301, 0.99, 2);

% plot the curve of convergence
% figure;
% plot(costs);
% title('Curve of Convergence');
% xlabel('Epoch'); ylabel('Cost');

%% Signal Decision and BER Calculation
% For the unequalized signal
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(ExtractedSignal, OriginalSignal);
fprintf('\nThe signal error before equalization\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% For the equalized signal
% [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
% fprintf('\nThe signal error after equalization\n');
% fprintf('Bit number num: %d \n', BitErrorNum);
% fprintf('SER: %e\n', SymErrorRate);
% fprintf('BER: %e\n', BitErrorRate);
