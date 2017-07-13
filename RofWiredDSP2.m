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
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

PathToSampledDataFolder = '.\Sampled Data\rof_20170608\';
FileList = {'-4.txt';};
% VolterraSetup = repmat([123, 45, 0, 0.003], length(FileList), 1);
BER = zeros(length(FileList), 1);
for i = 1 : length(FileList)
	PathToSampledData = strcat(PathToSampledDataFolder, FileList{i});

	SampledSignal = importdata(PathToSampledData);
	SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

	% Draw Eye Diagram
	% ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
	% step(ed1, SampledSignal);

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
	[EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalSignal, 'lms', 21, 0.01, 5);
	% figure;
	% plot(costs);
	% title('Curve of Convergence');
	% xlabel('Epoch'); ylabel('Cost');

	[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalSignal);
	fprintf('Equalization Setup: Linear LMS Channel Length is %d, alpha is 0.01\n', 21);
	fprintf('Bit number num: %d \n', BitErrorNum);
	fprintf('SER: %e\n', SymErrorRate);
	fprintf('BER: %e\n', BitErrorRate);
	toc
	BER(i) = BitErrorRate;
end

% x = [-3; -4; -5; -6; -7; -8; -9; -10; -11];
% plot(x, log10(BER))
% title('BER curve');
% xlabel('ROP/dBm'); ylabel('log10(BER)');
