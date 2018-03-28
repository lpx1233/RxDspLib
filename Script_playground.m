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
TxData = zeros(101, 1);
TxData(51) = 1;
TxSignal = reshape(repmat(TxData, 1, OverSamplingRatio)', [], 1);
TxSignal = lowPass10G(TxSignal);
TxSymNum = 1:0.25:101.75;
figure;plot(TxSymNum, TxSignal)
plotInFreq(TxSignal, SampleRate)

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
  % % % 20-order 10G Kaiser LPF
  % % % All frequency values are in GHz.
  % % Fs = 100;  % Sampling Frequency
  % % N    = 20;       % Order
  % % Fc   = 10;       % Cutoff Frequency
  % % flag = 'scale';  % Sampling Flag
  % % Beta = 1;        % Window Parameter
  % Hd = dsp.FIRFilter( ...
  %   'Numerator', [-6.64318883454626e-18 -0.01856525816569 ...
  %   -0.0352121682532305 -0.041700902608702 -0.0309951315585328 ...
  %   7.94840490930163e-18 0.0486756433241017 0.106687218997798 ...
  %   0.161838954883738 0.201391397662279 0.215760491436478 0.201391397662279 ...
  %   0.161838954883738 0.106687218997798 0.0486756433241017 ...
  %   7.94840490930163e-18 -0.0309951315585328 -0.041700902608702 ...
  %   -0.0352121682532305 -0.01856525816569 -6.64318883454626e-18]);

  % The following code was used to design the filter coefficients:
  %
  % N    = 20;   % Order
  % F3dB = 10;   % 3-dB Frequency
  % Fs   = 100;  % Sampling Frequency
  Hd = dsp.BiquadFilter( ...
      'Structure', 'Direct form II', ...
      'SOSMatrix', [1 2 1 1 -1.54670446522043 0.911831860114908; 1 2 1 1 ...
      -1.42280301693279 0.758681247520709; 1 2 1 1 -1.32091343081943 ...
      0.632738792885277; 1 2 1 1 -1.23786474339819 0.530084969790494; 1 2 1 1 ...
      -1.1710153071115 0.447454522282603; 1 2 1 1 -1.11823348180364 ...
      0.382212598225574; 1 2 1 1 -1.07784909638848 0.33229475262289; 1 2 1 1 ...
      -1.04859957636261 0.29614035756167; 1 2 1 1 -1.0295819074045 ...
      0.272633225955856; 1 2 1 1 -1.02021514617867 0.261055272351718], ...
      'ScaleValues', [0.09128184872362; 0.0839695576469805; ...
      0.0779563405164626; 0.0730550565980766; 0.0691098037927765; ...
      0.0659947791054833; 0.0636114140586013; 0.0618851952997645; ...
      0.060762829637839; 0.0602100315432626; 1]);
  y = step(Hd, double(x));
end
