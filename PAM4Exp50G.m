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
% eyediagram(SampledSignal(1:100000), 4*OverSamplingRatio, 2*OverSamplingRatio, 0.5*OverSamplingRatio);
% grid on;

%% Signal Synchronization and Extraction
[ExtractedSignal, OriginalSignal] = syncAndExtractSignal(SampledSignal, OriginalData, OverSamplingRatio);

[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(ExtractedSignal, OriginalSignal);
fprintf('The signal error before equalization\n');
fprintf('Bit number num: %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% ChannelLength = 35 : 2 : 41;
% BitErrorNum = zeros(length(ChannelLength), 1);
% for i = 1 : length(ChannelLength)
	tic
	%% Equalization Setup 1
	[EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 10, 123, 0.003, 43, [], 0, [], false, 0.001);
	% [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 123, 0.003, 50);
	figure;
	plot(costs);
	title('Curve of Convergence');
	xlabel('Epoch'); ylabel('Cost');

	% [EqualizedSignal, ChnlCoeffs, costs] = mlseEqualize(EqualizedSignal, OriginalSignal, 8);
	% figure;
	% plot(costs);
	% title('Curve of Convergence of MLSE for setup 1');
	% xlabel('Epoch'); ylabel('Cost');
	% Signal Decision and BER Calculation
	[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
	fprintf('\nEqualization Setup: Volterra LMS 1st order length = %d, 2nd order = %d, alpha = %f\n', 123, 43, 0.003);
	fprintf('Bit number num: %d \n', BitErrorNum);
	fprintf('SER: %e\n', SymErrorRate);
	fprintf('BER: %e\n', BitErrorRate);
	toc
% end
%
% figure;
% plot(ChannelLength, BitErrorNum);
% title('Curve of BitErrorNum');
% xlabel('FFE Channel Length'); ylabel('BitErrorNum');
