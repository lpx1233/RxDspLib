clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

OSCRate = 80e9;
DataRate = 28e9;
SampleRate = lcm(OSCRate, DataRate);
OverSamplingRatio = SampleRate / DataRate;
FileDir = '.\Sampled Data\50G PAM4\201709\56gpam4_10dbm\obtb\';
% ROP = -22 : -12;
% ed = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(SampledSignal), max(SampledSignal)]);
% step(ed, SampledSignal);
% BER = zeros(length(ROP), 5);
% BitError = zeros(length(ROP), 5);
% i = 1;
% for i = 1 : length(ROP)
    FileName = ['C2-12dBm00000.dat'];
    % FileName = [num2str(ROP(i)), 'dBm0.csv'];
    % FileName = '-10.txt';
    SampledSignal = importdata([FileDir, FileName]);
    SampledSignal = resample(SampledSignal, SampleRate, OSCRate);
    SampledSignal = (SampledSignal - mean(SampledSignal)) / std(SampledSignal);
    % % Sequence Extraction
    tic
    OriginalSignal = importdata('.\Original Data\Original_Data.txt');
    OriginalSignal = (OriginalSignal - 0.5) * 2;
    OriginalData_port1 = OriginalSignal;
    shiftnum = 19732;
    OriginalData_port2 = [-(OriginalSignal(end - shiftnum + 1 : end));
                          -(OriginalSignal(1 : end - shiftnum))];
    OriginalData = 2 * OriginalData_port1 + OriginalData_port2;
    % OriginalData = OriginalData_port1;
    OriginalDataUS = upsample(OriginalData, OverSamplingRatio);
    CorrelationResult = conv(SampledSignal, OriginalDataUS(end:-1:1), 'valid');
    figure;
    plot(CorrelationResult);
    title(FileName);
    toc
    [a, index] = max(CorrelationResult);
    ExtractedSignal = SampledSignal(index : OverSamplingRatio : index + length(OriginalData) * OverSamplingRatio - 1);

    % BER counting
    [BitErrorRate, SymErrorRate, BitErrorNum, OutputSignal] =  decisionAndCalcBerPAM4(ExtractedSignal, OriginalData);
    fprintf('The signal error before equalization\n');
    fprintf('Bit error num: %d\n', BitErrorNum);
    fprintf('SER: %e\n', SymErrorRate);
    fprintf('BER: %e\n', BitErrorRate);
    % BER(i, 5) = BitErrorRate;
    % BitError(i, 5) = BitErrorNum;

    % EQ and BER counting
    tic
    ChannelLen = 87;
    alpha = 0.0067;
    % for j = 1 : length(ChannelLen)
        [EqualizedSignal, w, costs] = linearFFEqualize(ExtractedSignal, OriginalData, 'lms', ChannelLen, alpha, 10);
        % figure;
        % plot(costs);
        % title('Curve of Convergence');
        % xlabel('Epoch'); ylabel('Cost');
        [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
        fprintf('Equalization Setup: Linear LMS Channel Length is %d, alpha is %f\n', ChannelLen, alpha);
        fprintf('Bit error num: %d\n', BitErrorNum);
        fprintf('SER: %e\n', SymErrorRate);
        fprintf('BER: %e\n', BitErrorRate);
        % BER(i, j) = BitErrorRate;
        % BitError(i, j) = BitErrorNum;
    % end
    toc

    % tic
    % ChanLen1st = 121;
    % ChanLen2nd = 21;
    % ChanLen3rd = 5;
    % alpha = 0.003;
    % [EqualizedSignal, w, costs] = volterraFFEqualize(ExtractedSignal, OriginalData, 'lms', 10, ChanLen1st, alpha, ChanLen2nd, [], ChanLen3rd);
    % figure;
    % plot(costs);
    % title('Curve of Convergence');
    % xlabel('Epoch'); ylabel('Cost');
    %
    % EqualizedSignalUS = resample(EqualizedSignal, OverSamplingRatio, 1);
    % ed1 = comm.EyeDiagram('DisplayMode','2D color histogram','OversamplingMethod','Input interpolation', 'SamplesPerSymbol', OverSamplingRatio, 'YLimits', [min(EqualizedSignalUS), max(EqualizedSignalUS)]);
    % step(ed1, EqualizedSignalUS);
    %
    % [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(EqualizedSignal, OriginalData);
    % fprintf('Equalization Setup: Volterra LMS Channel Length Setup is [%d %d %d], alpha is %f\n', ChanLen1st, ChanLen2nd, ChanLen3rd, alpha);
    % fprintf('Bit error num: %d\n', BitErrorNum);
    % fprintf('SER: %e\n', SymErrorRate);
    % fprintf('BER: %e\n', BitErrorRate);
    % toc
% end
