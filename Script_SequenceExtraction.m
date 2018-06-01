clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

OSCRate = 80e9;
DataRate = 25e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;
ProjectDir = '.\Sampled Data\201805lt\';
FileDir = [ProjectDir, '20180519\25G PAM16\soa\btb\prbs15\'];
ROP = -7;
i = 0;
% for ROP = -7 : -2 : -15
  for i = 0 : 9
    FileName = [num2str(ROP), 'dBm', num2str(i), '.txt'];
    FilePath = [FileDir, 'raw\', FileName];
    fprintf(['Processing ', replace(FilePath, '\', '\\'), ' ...\n']);
    SampledSignal = importdata(FilePath);
    SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);
    % plotInFreq(SampledSignal, OSCRate);
    SampledSignal = LPF(SampledSignal);
    % plotInFreq(SampledSignal, OSCRate);
    SampledSignal = resample(SampledSignal, SampleRate, OSCRate);

    % OriginalData = importdata([ProjectDir, 'pam16_', num2str(i), '.csv']);
    OriginalData = importdata([ProjectDir, 'prbs15_pam16.csv']);
    OriginalData = (OriginalData - mean(OriginalData)) / std(OriginalData);

    tic
    OriginalDataUS = upsample(OriginalData, OverSamplingRatio);
    CorrelationResult = conv(SampledSignal, OriginalDataUS(end:-1:1), 'valid');
    figure;
    plot(CorrelationResult);
    title(FileName);
    toc
    [a, index] = max(CorrelationResult);

    ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);
    csvwrite([FileDir, 'filtered\', num2str(ROP), 'dBm', num2str(i), '.csv'], ExtractedSignal);
  end
% end

%% Aux Functions
function y = plotInFreq(in, Fs)
  % plotInFreq: plot signal in frequency domain
  L = length(in);
  Y = fft(in);
  M2 = abs(Y / L);
  M1 = M2(1 : L/2 + 1);
  M1(2 : end - 1) = 2 * M1(2 : end - 1);
  P2 = angle(Y / L);
  P1 = P2(1 : L/2 + 1) / pi;
  f = Fs * (0 : (L / 2)) / L;
  figure;
  % subplot(2, 1, 1);
  plot(f, 20 * log10(M1))
  % subplot(2, 1, 2);
  % plot(f, P1)
end

function y = LPF(x)
  % All frequency values are in GHz.
  Fs = 80;  % Sampling Frequency
  Fpass = 18;              % Passband Frequency
  Fstop = 22;              % Stopband Frequency
  Dpass = 0.057501127785;  % Passband Ripple
  Dstop = 0.0001;          % Stopband Attenuation
  dens  = 20;              % Density Factor
  % Calculate the order from the parameters using FIRPMORD.
  [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
  % Calculate the coefficients using the FIRPM function.
  b  = firpm(N, Fo, Ao, W, {dens});
  Hd = dsp.FIRFilter('Numerator', b);
  y = step(Hd,double(x));
end
