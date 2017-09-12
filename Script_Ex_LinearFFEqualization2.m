%% This example shows how to equalize all the dat file in the same directory
%% and then draw a ber plot for the dat files with 1 points a file.

% TODO The ber plot is not done yet

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
SampleDir = '.\Sampled Data\40km+FBG+FILTER+0\';
files = dir(SampleDir);
for i = 1 : length(files) - 2
	SampledSignal = importdata(strcat(SampleDir,files(i + 2).name));
	% eyediagram(SampledSignal(1:100000), 4*OverSamplingRatio, 2*OverSamplingRatio, 0.5*OverSamplingRatio);
	% grid on;

	%% Signal Synchronization and Extraction
	[ExtractedSignal(:,i), OriginalSignal(:,i)] = syncAndExtractSignal(SampledSignal, OriginalData, OverSamplingRatio);

	%% LMS Equalization
	% 101-tap FFE and training for 20 epochs
	[equalizedSignal(:,i), w] = linearFFEqualize(ExtractedSignal(:,i), OriginalSignal(:,i), 'lms', 101, 0.01, 20);
	% Equalization Result Visulization
	% eyediagram(equalizedSignal, 4*OverSamplingRatio, 2*OverSamplingRatio, 0.5*OverSamplingRatio);
	% grid on;

	%% Signal Decision and BER Calculation
	[BitErrorRate(i), SymErrorRate(i), BitErrorNum(i)] = decisionAndCalcBerPAM4(equalizedSignal(:,i), (OriginalSignal(:,i) + 3) / 2);
	fprintf('The %d th Bit eror numbers: %d \n', i, BitErrorNum(i));
	fprintf('SER: %e\n', SymErrorRate(i));
	fprintf('BER: %e\n\n', BitErrorRate(i));
end
	plot(log10(BitErrorRate))