clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

DataLength = 500000;
SampleRate = 100e9;
SymbolRate = 25e9;
OverSamplingRatio = SampleRate / SymbolRate;
rng('shuffle');
% Tx signal generation
PAM4Data = randi(4, DataLength, 1) - 1;
TxData = [PAM4Data(end - 99: end); PAM4Data; PAM4Data(1: 100)];
% TxSignal = resample(TxData, SampleRate, SymbolRate);
TxSignal = reshape(repmat(TxData, 1, OverSamplingRatio)', [], 1);
TxSignal = (TxSignal - mean(TxSignal)) / std(TxSignal);
plotInFreq(TxSignal, SampleRate)

% Transmit through a band-limited AWGN channel
RxSignal = TxSignal;
RxSignal = zeros(size(TxSignal));
% RxSignal = awgn(RxSignal, 22.5);
% RxSignal = lowPass10G(RxSignal);
RxSignal = awgn(RxSignal, 20);
% RxSignal = lowPass30G(RxSignal);
% RxSignal = (RxSignal - mean(RxSignal)) / std(RxSignal);
% RxSignal = lowPass10G(RxSignal);
plotInFreq(RxSignal, SampleRate)

% % TODO remove this
% ed0 = comm.EyeDiagram('DisplayMode', '2D color histogram', 'OversamplingMethod', 'Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(RxSignal), max(RxSignal)]);
% step(ed0, RxSignal);

% TODO Matched Filtering
RxSignal = filter(ones(OverSamplingRatio, 1), 1, RxSignal);
plotInFreq(RxSignal, SampleRate)

% % Extracted transmitted signal
% tic
% OriginalDataUS = upsample(PAM4Data, OverSamplingRatio);
% CorrelationResult = conv(RxSignal, OriginalDataUS(end:-1:1), 'valid');
% figure;
% plot(CorrelationResult);
% toc
% [a, index] = max(CorrelationResult);
%
% ExtractedSignal = RxSignal(index : index + DataLength * OverSamplingRatio - 1);
% ed = comm.EyeDiagram('DisplayMode', '2D color histogram', 'OversamplingMethod', 'Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(ExtractedSignal), max(ExtractedSignal)]);
% step(ed, ExtractedSignal);
%
% ExtractedSignal = RxSignal(index : OverSamplingRatio : index + DataLength * OverSamplingRatio - 1);
% Dataset = zeros(DataLength - 100, 101);
% Noise = awgn(zeros(size(ExtractedSignal)), 25);
% NoiseDataset = zeros(DataLength - 100, 101);;
% for i = 1 : length(Dataset)
%   Dataset(i, :) = ExtractedSignal(i : i + 100, :)';
%   NoiseDataset(i, :) = Noise(i : i + 100, :)';
% end
% figure;
% subplot(2, 1, 1); plot(svd(Dataset))
% subplot(2, 1, 2); plot(svd(NoiseDataset))

% x = zeros(length(ExtractedSignal) - 101 + 1, 101);
% for i = 1 : length(ExtractedSignal) - 101 + 1
%   x(i, :) = ExtractedSignal(i : i + 101 - 1)';
% end
% Sigma = (1 / size(x, 1)) * x' * x;
% [U, S, V] = svd(Sigma);
% k = 70;
% z = x * U(:, 1 : k);
% x_approx = z * U(:, 1 : k)';
% figure;plot(sum(S));
% sum(sum((x - x_approx) .^ 2)) / sum(sum((x) .^ 2))
% csvwrite(['.\Sampled Data\Simulation\', 'rx.csv'], ExtractedSignal);
% csvwrite(['.\Sampled Data\Simulation\', 'tx.csv'], PAM4Data);

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
  % % 20-order 10G Kaiser LPF
  % % All frequency values are in GHz.
  % Fs = 100;  % Sampling Frequency
  % N    = 20;       % Order
  % Fc   = 10;       % Cutoff Frequency
  % flag = 'scale';  % Sampling Flag
  % Beta = 1;        % Window Parameter
  Hd = dsp.FIRFilter( ...
    'Numerator', [-6.64318883454626e-18 -0.01856525816569 ...
    -0.0352121682532305 -0.041700902608702 -0.0309951315585328 ...
    7.94840490930163e-18 0.0486756433241017 0.106687218997798 ...
    0.161838954883738 0.201391397662279 0.215760491436478 0.201391397662279 ...
    0.161838954883738 0.106687218997798 0.0486756433241017 ...
    7.94840490930163e-18 -0.0309951315585328 -0.041700902608702 ...
    -0.0352121682532305 -0.01856525816569 -6.64318883454626e-18]);
  y = step(Hd, double(x));

  % % A butterworth LPF simulating DML/PD inband frequency fading
  % Fpass = 1;    % Passband Frequency
  % Fstop = 20;   % Stopband Frequency
  % Apass = 1;    % Passband Ripple (dB)
  % Astop = 10;   % Stopband Attenuation (dB)
  % Fs    = 100;  % Sampling Frequency
  % h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);
  % Hd = design(h, 'butter', ...
  %     'MatchExactly', 'stopband', ...
  %     'SystemObject', true);
  % x = step(Hd, double(x));
  % % 4th-order Bessel LPF
  % % All frequency values are in GHz.
  % Nb   = 4;    % Numerator Order
  % Na   = 4;    % Denominator Order
  % F3dB = 12;   % 3-dB Frequency
  % Fs   = 100;  % Sampling Frequency
  % h = fdesign.lowpass('nb,na,f3db', Nb, Na, F3dB, Fs);
  % Hd = design(h, 'butter', 'SystemObject', true);
  % y = step(Hd, double(x));
end

function y = lowPass25G(x)
  % % All frequency values are in GHz.
  % Fs = 100;  % Sampling Frequency
  % N    = 20;       % Order
  % Fc   = 25;       % Cutoff Frequency
  % flag = 'scale';  % Sampling Flag
  % Beta = 1;        % Window Parameter
  Hd = dsp.FIRFilter( ...
    'Numerator', [1.50218896096586e-17 0.0285686861166202 ...
    -1.64045849115145e-17 -0.0396594867673816 1.75232583993696e-17 ...
    0.0587052701799814 -1.83460182902814e-17 -0.101464478833822 ...
    1.88493140942653e-17 0.309906147011409 0.487887724586386 ...
    0.309906147011409 1.88493140942653e-17 -0.101464478833822 ...
    -1.83460182902814e-17 0.0587052701799814 1.75232583993696e-17 ...
    -0.0396594867673816 -1.64045849115145e-17 0.0285686861166202 ...
    1.50218896096586e-17]);
  y = step(Hd,double(x));
end

function y = lowPass5G(x)
  % % All frequency values are in GHz.
  % Fs = 100;  % Sampling Frequency
  % N    = 20;       % Order
  % Fc   = 5;        % Cutoff Frequency
  % flag = 'scale';  % Sampling Flag
  % Beta = 1;        % Window Parameter
  Hd = dsp.FIRFilter( ...
      'Numerator', [2.71637240797039e-18 0.00798192017345794 ...
      0.017797042461414 0.0290094491915146 0.0410132167097605 ...
      0.0530776687430912 0.0644083314945248 0.0742175172547654 ...
      0.0817971427167792 0.0865859254643814 0.0882235715806219 ...
      0.0865859254643814 0.0817971427167792 0.0742175172547654 ...
      0.0644083314945248 0.0530776687430912 0.0410132167097605 ...
      0.0290094491915146 0.017797042461414 0.00798192017345794 ...
      2.71637240797039e-18]);
  y = step(Hd,double(x));
end

function y = lowPass30G(x)
  Fpass = 30;   % Passband Frequency
  Fstop = 33;   % Stopband Frequency
  Apass = 1;    % Passband Ripple (dB)
  Astop = 80;   % Stopband Attenuation (dB)
  Fs    = 100;  % Sampling Frequency
  h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);
  Hd = design(h, 'butter', ...
      'MatchExactly', 'stopband', ...
      'SystemObject', true);
  y = step(Hd, double(x));
end
