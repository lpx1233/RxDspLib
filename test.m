clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

RawSignal = importdata('.\Sampled Data\RoF\wireless\rof_btb_eml20170713.dat');

SampleRate = 600e9;
OSCRate = 120e9;
DataRate = 12.5e9;
OverSamplingRatio = SampleRate / DataRate;

SampledSignal = 100 * (RawSignal(:, 2) - mean(RawSignal(:, 2)));
SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

t = RawSignal(:, 1);
t = interp1(1:length(t), t, 1:OSCRate/SampleRate:length(t)+1);
t = t';
t = t(1:end-1);

index = find(abs(t) < 1e-9);
lo = sin(2*pi*25.5e9*t);

plot(t(index), SampledSignal(index), t(index), lo(index), 'r')
