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
SampledSignal = importdata('.\Sampled Data\10dml25km10Apd\-16dsf.txt');
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
eyediagram(SampledSignal(1:100000), 4*OverSamplingRatio, 2*OverSamplingRatio, 0.5*OverSamplingRatio);
grid on;

%% Signal Synchronization and Extraction
[ExtractedSignal, OriginalSignal] = syncAndExtractSignal(SampledSignal, OriginalData, OverSamplingRatio);

[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(ExtractedSignal, OriginalSignal);
fprintf('The signal error before equalization\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

%% Equalization Setup 1
[EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalSignal, 'rls', 3, 15, 0.999, [], [], [], [], false, 0.001);
figure;
plot(costs);
title('Curve of Convergence of Volterra for setup 1');
xlabel('Epoch'); ylabel('Cost');

[EqualizedSignal, ChnlCoeffs, costs] = mlseEqualize(EqualizedSignal, OriginalSignal, 8);
figure;
plot(costs);
title('Curve of Convergence of MLSE for setup 1');
xlabel('Epoch'); ylabel('Cost');
% Signal Decision and BER Calculation
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
fprintf('\nEqualization Setup 1: RLS Volterra length = 15, lambda = 0.999; MLSE states = 8\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

%% Equalization Setup 2
[EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 10, 15, 0.03, [], [], [], [], false);
figure;
plot(costs);
title('Curve of Convergence of Volterra for setup 2');
xlabel('Epoch'); ylabel('Cost');

[EqualizedSignal, ChnlCoeffs, costs] = mlseEqualize(EqualizedSignal, OriginalSignal, 8);
figure;
plot(costs);
title('Curve of Convergence of MLSE for setup 2');
xlabel('Epoch'); ylabel('Cost');
% Signal Decision and BER Calculation
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
fprintf('\nEqualization Setup 2: LMS Volterra length = 15, alpha = 0.03; MLSE states = 8\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

%% Equalization Setup 3
[EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalSignal, 'rls', 3, 17, 0.999, [], [], [], [], false, 0.001);
figure;
plot(costs);
title('Curve of Convergence of Volterra for setup 3');
xlabel('Epoch'); ylabel('Cost');

[EqualizedSignal, ChnlCoeffs, costs] = mlseEqualize(EqualizedSignal, OriginalSignal, 8);
figure;
plot(costs);
title('Curve of Convergence of MLSE for setup 3');
xlabel('Epoch'); ylabel('Cost');
% Signal Decision and BER Calculation
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
fprintf('\nEqualization Setup 3: RLS Volterra length = 17, lambda = 0.999; MLSE states = 8\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

%% Equalization Setup 4
[EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 10, 17, 0.03, [], [], [], [], false);
figure;
plot(costs);
title('Curve of Convergence of Volterra for setup 4');
xlabel('Epoch'); ylabel('Cost');

[EqualizedSignal, ChnlCoeffs, costs] = mlseEqualize(EqualizedSignal, OriginalSignal, 8);
figure;
plot(costs);
title('Curve of Convergence of MLSE for setup 4');
xlabel('Epoch'); ylabel('Cost');
% Signal Decision and BER Calculation
[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
fprintf('\nEqualization Setup 2: LMS Volterra length = 17, alpha = 0.03; MLSE states = 8\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);
