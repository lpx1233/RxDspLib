clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

OSCRate = 80e9;
DataRate = 25e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;
FileDir = '.\Sampled Data\201801\O-band DML\40G PD\obtb\';
% FileDir = '.\Sampled Data\50G PAM4\201709\56gpam4_10dbm\25km\';
% ROP = -22:-13;
ROP = 5;
i = 0;
% FileDirList = {'.\Sampled Data\201710\112G\56G_PAM4\25G_APD\50km+1700\', '.\Sampled Data\201710\112G\56G_PAM4\25G_APD\50km+850\', '.\Sampled Data\201710\112G\56G_PAM4\25G_APD\50km\', '.\Sampled Data\201710\112G\56G_PAM4\20G_PIN_EDFA\50km\', '.\Sampled Data\201710\112G\56G_PAM4\20G_PIN_EDFA\50km+850\', '.\Sampled Data\201710\112G\56G_PAM4\20G_PIN_EDFA\50km+1700\', '.\Sampled Data\201710\112G\56G_PAM4\20G_PIN_EDFA\obtb\'};
% for k = 1 : 7
%   FileDir = FileDirList{k};
  % for ROP = -3 : 2 : 3
  %   for i = 0 : 4
      FileName = ['C2', num2str(ROP), 'dBm0000', num2str(i), '.dat'];
      % FileName = ['C2obtb0000', num2str(i), '.dat'];
      fprintf(['Processing ', replace([FileDir, 'raw\'], '\', '\\'), FileName, ' ...\n']);
      % FileName = 'C2ebtb00000.dat';
      SampledSignal = importdata([FileDir, 'raw\', FileName]);
      SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
      SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);

      OriginalData = importdata(['.\Sampled Data\201801\pam4_', num2str(i), '.csv']);
      % OriginalData = importdata('.\Sampled Data\50G PAM4\201709\56gpam4_10dbm\25km\extracted\prbs15pam4.csv');
      % OriginalData = OriginalData';
      % OriginalData = 3 - OriginalData;
      % csvwrite([FileDir, 'extracted\pam4_', num2str(i), '.csv'], OriginalData);
      OriginalData = (OriginalData - mean(OriginalData)) / std(OriginalData);
      % OriginalSignal = (OriginalSignal - 0.5) * 2;
      % OriginalData_port1 = OriginalSignal;
      % shiftnum = 13035;
      % OriginalData_port2 = [-(OriginalSignal(shiftnum + 1 : end));
      % 				-(OriginalSignal(1 : shiftnum))];
      % OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
      % tic
      % CorrelationResult = zeros(length(SampledSignal) - OverSamplingRatio * length(OriginalData) + 1, 1);
      % parfor i = 1 : length(CorrelationResult)
      %   CorrelationResult(i) = sum(SampledSignal(i : OverSamplingRatio : i + OverSamplingRatio * length(OriginalData) - 1) .* OriginalData);
      % end
      % figure;
      % plot(CorrelationResult);
      % title(FileName);
      % toc

      % plotInFreq(SampledSignal, SampleRate);
      % SampledSignal = lowPass10G(SampledSignal);
      % SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);
      plotInFreq(SampledSignal, SampleRate);

      % TODO remove this
      ed0 = comm.EyeDiagram('DisplayMode', '2D color histogram', 'OversamplingMethod', 'Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
      step(ed0, SampledSignal);
      % TODO Matched Filtering
      SampledSignal = filter(ones(OverSamplingRatio, 1), 1, SampledSignal);
      SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);
      plotInFreq(SampledSignal, SampleRate);

      tic
      OriginalDataUS = upsample(OriginalData, OverSamplingRatio);
      CorrelationResult = conv(SampledSignal, OriginalDataUS(end:-1:1), 'valid');
      figure;
      plot(CorrelationResult);
      title(FileName);
      toc
      [a, index] = max(CorrelationResult);

      % ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);
      % csvwrite([FileDir, 'extracted\', num2str(ROP), 'dBm', num2str(i), '.csv'], ExtractedSignal);

      ExtractedSignal = SampledSignal(index : index + length(OriginalData) * OverSamplingRatio - 1);
      ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Histogram interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(ExtractedSignal), max(ExtractedSignal)]);
      step(ed, ExtractedSignal);
  %   end
  % end
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

function y = lowPass10G(x)
  % All frequency values are in GHz.
  Fs = 560;  % Sampling Frequency
  Fpass = 10;              % Passband Frequency
  Fstop = 15;              % Stopband Frequency
  Dpass = 0.057501127785;  % Passband Ripple
  Dstop = 0.01;            % Stopband Attenuation
  flag  = 'scale';         % Sampling Flag
  % Calculate the order from the parameters using KAISERORD.
  [N, Wn, BETA, TYPE] = kaiserord([Fpass Fstop]/(Fs/2), [1 0], [Dstop Dpass]);
  % Calculate the coefficients using the FIR1 function.
  b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
  Hd = dsp.FIRFilter('Numerator', b);
  y = step(Hd,double(x));
end
