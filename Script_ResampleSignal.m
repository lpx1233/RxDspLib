clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

OSCRate = 80e9;
DataRate = 28e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;
FileDir = '.\Sampled Data\201710\28G_PAM4\obtb\';
files = dir([FileDir, 'extracted\']);

for i = 1 : length(files) - 2
  SampledSignal = importdata([FileDir, 'extracted\', files(i).name]);
  ExtractedSignal = SampledSignal(1 : OverSamplingRatio : end);
  csvwrite([FileDir, 'extracted_1sps\', files(i).name], ExtractedSignal);
end

load splat
sound(y,Fs)
