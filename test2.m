clear all;
close all;
clc;
%% Change dir and recovery workspace
cd(fileparts(which(mfilename)));
load('.\Sampled Data\RoF\wireless\40km\20170727tdc\rof_3dBm_corr.mat')

[a, index] = max(CorrelationResult);

ExtractedSignal = r(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);

% threshold_lo = 0:1/40:28/80;
% threshold_mi = 3/8:1/40:5/8;
% threshold_hi = 52/80:1/40:1;
% BER = zeros(length(threshold_lo), length(threshold_mi), length(threshold_hi));
% BitError = BER;
% for i = 1 : length(threshold_lo)
%   for j = 1 : length(threshold_mi)
%     for k = 1 : length(threshold_hi)
%       threshold = [threshold_lo(i), threshold_mi(j), threshold_hi(k)];
%       [BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData, threshold);
%       fprintf('The signal error before equalization with threshold [%f, %f, %f]\n', threshold);
%       fprintf('Bit error num: %d\n', BitErrorNum);
%       fprintf('SER: %e\n', SymErrorRate);
%       fprintf('BER: %e\n', BitErrorRate);
%       BitError(i, j, k) = BitErrorNum;
%       BER(i, j, k) = BitErrorRate;
%     end
%   end
% end

[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
fprintf('The signal error before equalization\n');
fprintf('Bit error num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);


% % linear FFE
%
% %% TODO remove this if not using T/2-spaced FFE
% ExtractedSignal = r(index : OverSamplingRatio/2 : index + length(OriginalData) * OverSamplingRatio - 1);
% OriginalDataUS = interp1(1:length(OriginalData), OriginalData, 1 : 0.5 : length(OriginalData)+0.5);
% OriginalDataUS = OriginalDataUS';
% OriginalDataUS(end) = 0;
%
% ChannelLen = 101:20:1001;
% alpha = [0.001; 0.003];
% BER = zeros(length(ChannelLen), length(alpha));
% BitError = zeros(length(ChannelLen), length(alpha));
% for i = 1 : length(ChannelLen)
%   for j = 1 : length(alpha)
%     [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalDataUS, 'lms', ChannelLen(i), alpha(j), 5);
%     % figure;
%     % plot(costs);
%     % title('Curve of Convergence');
%     % xlabel('Epoch'); ylabel('Cost');
%
%     %% TODO remove this if not using T/2-spaced FFE
%     EqualizedSignal = EqualizedSignal(1:2:end);
%
%     [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
%     fprintf('Equalization Setup: Linear LMS Channel Length is %d, alpha is %f\n', ChannelLen(i), alpha(j));
%     fprintf('Bit error num: %d\n', BitErrorNum);
%     fprintf('SER: %e\n', SymErrorRate);
%     fprintf('BER: %e\n', BitErrorRate);
%     BitError(i, j) = BitErrorNum;
%     BER(i, j) = BitErrorRate;
%   end
% end
% meshz(log10(alpha), ChannelLen, -log10(BER));
% title('BER vs ChannelLen & alpha')
% xlabel('log10(alpha)') % x-axis label
% ylabel('ChannelLen') % y-axis label
% zlabel('-log10(BER)') % y-axis label

% Volterra
% % ChanLen1st = 101:10:301;
% % ChanLen2nd = [0, 1:2:21];
% % ChanLen3rd = [0, 1:2:7];
% % alpha = [0.01; 0.003; 0.001];
% % BER = zeros(length(ChanLen1st), length(ChanLen2nd), length(ChanLen3rd), length(alpha));
% % for i = 1 : length(ChanLen1st)
% %   for j = 1 : length(ChanLen2nd)
% %     for k = 1 : length(ChanLen3rd)
% %       for l = 1 : length(alpha)
          tic

          % TODO remove these parameter
          ChanLen1st = 301;
          ChanLen2nd = 33;
          ChanLen3rd = 11;
          alpha = 0.003;

          [EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalData, 'lms', 10, ChanLen1st, alpha, ChanLen2nd, [], ChanLen3rd);
          figure;
          plot(costs);
          title('Curve of Convergence');
          xlabel('Epoch'); ylabel('Cost');

          EqualizedSignalUS = resample(EqualizedSignal, 48, 1);
          ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', 48, 'YLimits', [min(EqualizedSignalUS), max(EqualizedSignalUS)]);
          step(ed, EqualizedSignalUS);

          [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
          fprintf('Equalization Setup: Volterra LMS Channel Length Setup is [%d %d %d], alpha is %f\n', ChanLen1st, ChanLen2nd, ChanLen3rd, alpha);
          fprintf('Bit error num: %d\n', BitErrorNum);
          fprintf('SER: %e\n', SymErrorRate);
          fprintf('BER: %e\n', BitErrorRate);
          % BER(i, j, k, l) = BitErrorRate;
          toc
% %       end
% %     end
% %   end
% % end
% % meshz(log10(alpha), ChannelLen, -log10(BER));
% % title('BER vs ChannelLen & alpha')
% % xlabel('log10(alpha)') % x-axis label
% % ylabel('ChannelLen') % y-axis label
% % zlabel('-log10(BER)') % y-axis label
