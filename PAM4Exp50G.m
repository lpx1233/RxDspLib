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

PathToSampledDataFolder = '.\Sampled Data\20170428\wSOA wDSF 25km\soa12\';
FileList = {'6_0_7_-10.txt'; '6_0_7_-13.txt'; '6_0_7_-15.txt'; '6_0_7_-17.txt'; '6_0_7_-19.txt';};
VolterraSetup = repmat([123, 45, 0, 0.003], length(FileList), 1);

for i = 1 : length(FileList)
	PathToSampledData = strcat(PathToSampledDataFolder, FileList{i});

	SampledSignal = importdata(PathToSampledData);
	SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

	%% Signal Synchronization and Extraction
	[ExtractedSignal, OriginalSignal] = syncAndExtractSignal(SampledSignal, OriginalData, OverSamplingRatio);

	[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] = decisionAndCalcBerPAM4(ExtractedSignal, OriginalSignal);
	if i ~= 1
		fprintf('\n');
	end
	fprintf('The signal error before equalization for %s\n', PathToSampledData);
	fprintf('Bit number num: %d \n', BitErrorNum);
	fprintf('SER: %e\n', SymErrorRate);
	fprintf('BER: %e\n', BitErrorRate);

	tic
	[EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 10, VolterraSetup(i, 1), VolterraSetup(i, 4), VolterraSetup(i, 2), [], VolterraSetup(i, 3), [], false);
	figure;
	plot(costs);
	title('Curve of Convergence');
	xlabel('Epoch'); ylabel('Cost');

	[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
	fprintf('Equalization Setup: Volterra LMS Channel Length Setup is [%d %d %d], alpha is %f\n', VolterraSetup(i, :));
	fprintf('Bit number num: %d \n', BitErrorNum);
	fprintf('SER: %e\n', SymErrorRate);
	fprintf('BER: %e\n', BitErrorRate);
	toc
end
