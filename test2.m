clear all;
close all;
clc;
%% Change dir and recovery workspace
cd(fileparts(which(mfilename)));
load('.\Sampled Data\RoF\wireless\rof_btb_eml_corr.mat')

[a, index] = max(CorrelationResult);

ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);

[BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
fprintf('\nThe signal error before equalization\n');
fprintf('Bit number num: %d\n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

% % linear FFE
% ChannelLen = 101:10:301;
% alpha = [0.01; 0.003; 0.001];
% BER = zeros(length(ChannelLen), length(alpha));
% for i = 1 : length(ChannelLen)
%   for j = 1 : length(alpha)
%     [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalData, 'lms', ChannelLen(i), alpha(j), 5);
%     % figure;
%     % plot(costs);
%     % title('Curve of Convergence');
%     % xlabel('Epoch'); ylabel('Cost');
%     [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
%     fprintf('Equalization Setup: Linear LMS Channel Length is %d, alpha is %f\n', ChannelLen(i), alpha(j));
%     fprintf('Bit number num: %d\n', BitErrorNum);
%     fprintf('SER: %e\n', SymErrorRate);
%     fprintf('BER: %e\n', BitErrorRate);
%     BER(i, j) = BitErrorRate;
%   end
% end
% meshz(log10(alpha), ChannelLen, -log10(BER));
% title('BER vs ChannelLen & alpha')
% xlabel('log10(alpha)') % x-axis label
% ylabel('ChannelLen') % y-axis label
% zlabel('-log10(BER)') % y-axis label

ChanLen1st = 101:10:301;
ChanLen2nd = [0, 1:2:21];
ChanLen3rd = [0, 1:2:7];
alpha = [0.01; 0.003; 0.001];
BER = zeros(length(ChanLen1st), length(ChanLen2nd), length(ChanLen3rd), length(alpha));
for i = 1 : length(ChanLen1st)
  for j = 1 : length(ChanLen2nd)
    for k = 1 : length(ChanLen3rd)
      for l = 1 : length(alpha)
        tic
        [EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalData, 'lms', 10, ChanLen1st, alpha, ChanLen2nd, [], ChanLen3rd);
        % figure;
        % plot(costs);
        % title('Curve of Convergence');
        % xlabel('Epoch'); ylabel('Cost');
        [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
        fprintf('Equalization Setup: Volterra LMS Channel Length Setup is [%d %d %d], alpha is %f\n', ChanLen1st(i), ChanLen2nd(j), ChanLen3rd(k), alpha(l));
        fprintf('Bit number num: %d\n', BitErrorNum);
        fprintf('SER: %e\n', SymErrorRate);
        fprintf('BER: %e\n', BitErrorRate);
        BER(i, j, k, l) = BitErrorRate;
        toc
      end
    end
  end
end
% meshz(log10(alpha), ChannelLen, -log10(BER));
% title('BER vs ChannelLen & alpha')
% xlabel('log10(alpha)') % x-axis label
% ylabel('ChannelLen') % y-axis label
% zlabel('-log10(BER)') % y-axis label
