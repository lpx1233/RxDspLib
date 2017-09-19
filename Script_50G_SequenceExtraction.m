clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

OSCRate = 80e9;
DataRate = 28e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;
FileDir = '.\Sampled Data\50G PAM4\201709\56gpam4_10dbm\25km\';
FileList = dir(FileDir);
ROP = '-19dBm';
i = 0;
for i = 0 : 9
  FileName = ['C2', ROP, '0000', num2str(i), '.dat'];
  SampledSignal = importdata([FileDir, FileName]);
  SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
  SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);

  % ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
  % step(ed, SampledSignal);

  tic
  OriginalSignal = importdata('.\Original Data\Original_Data.txt');
  % OriginalData = OriginalData';
  % OriginalData = (OriginalData - mean(OriginalData)) * 6;
  OriginalSignal = (OriginalSignal - 0.5) * 2;
  OriginalData_port1 = OriginalSignal;
  shiftnum = 13035;
  OriginalData_port2 = [-(OriginalSignal(shiftnum + 1 : end));
						-(OriginalSignal(1 : shiftnum))];
  OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
  CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
  parfor i = 1 : length(CorrelationResult)
    CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
  end
  figure;
  plot(CorrelationResult);
  title(FileName);
  toc

  [a, index] = max(CorrelationResult);
  ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);
  [BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
  fprintf('The signal error before equalization\n');
  fprintf('Bit error num: %d\n', BitErrorNum);
  fprintf('SER: %e\n', SymErrorRate);
  fprintf('BER: %e\n', BitErrorRate);

  ExtractedSignal = SampledSignal(index : index + length(OriginalData) * OverSamplingRatio - 1);
  csvwrite(['.\Sampled Data\50G PAM4\201709\56gpam4_10dbm\25km\extracted\', ROP, num2str(i), '.csv'], ExtractedSignal);
end
