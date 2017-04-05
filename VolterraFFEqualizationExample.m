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
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;
% importing and eyediagram drawing
SampledSignal = importdata('.\Sampled Data\40km+FBG+FILTER+-680\F2_00011.dat');
% eyediagram(SampledSignal(1:100000), 4*OverSamplingRatio, 2*OverSamplingRatio, 0.5*OverSamplingRatio);
% grid on;

%% Signal Synchronization and Extraction
[ExtractedSignal, OriginalSignal] = syncAndExtractSignal(SampledSignal, OriginalData, OverSamplingRatio);

%% Volterra FFE Equalization
[EqualizedSignal, w, costs] = volLMSFFEqualize(ExtractedSignal, OriginalSignal, 5, 0.01, 10);

% plot the curve of convergence
figure;
plot(costs);
title('Curve of Convergence');
xlabel('Epoch'); ylabel('Cost');

%% Signal Decision and BER Calculation
% For the unequalized signal
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(ExtractedSignal, OriginalSignal);
fprintf('\nThe signal error before equalization\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% For the equalized signal
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
fprintf('\nThe signal error after equalization\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);