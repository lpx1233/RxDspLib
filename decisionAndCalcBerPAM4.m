function [BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(InputSignal, OriginalData, threshold)
	% This function performs the PAM4 decision and bit error ratio calculation for PAM4.
	% First the PAM4 decision will be performed and then the error counting and error ratio
	%	calculation. The %threshold% for PAM4 decision is a optional parameter, which have the
	% default value [0.25; 0.5; 0.75], allowing outside to perform some optimization on
	% the decision threshold of PAM4.
	%
	% input:
	%     InputSignal
	%       The input signal to be decided and error-calculated.
	%     OriginalData
	%       The origin data to be compared to the decided input signal to get the errors.
	%     threshold (optional)
	%       The PAM4 decision threshold, which can be adjusted outside this function.
	%       Default: [0.25; 0.5; 0.75]
	% output:
	%     BitErrorRate
	%       The bit error ratio of the %InputSignal% comparing to %OriginalData%.
	%     SymErrorRate
	%       The symbol error ratio of the %InputSignal% comparing to %OriginalData%.
	%     BitErrorNum
	%       The number of bit error, which is 1/sym when only 1 bit changes and 2/sym 
	%       when both bits change.
	
	%% Parameters Checking
	narginchk(2, 3);
	
	if ~exist('threshold','var') || isempty(threshold)
		threshold = [0.25; 0.5; 0.75];
	end
	
	%% Input Signal Normalization
	InputSignal = InputSignal - min(InputSignal);
	InputSignal = InputSignal / max(InputSignal);
	OriginalData = OriginalData - min(OriginalData);
	OriginalData = (OriginalData / max(OriginalData)) * 3;
	
	%% Input Signal Decision
	InputSignal(find(InputSignal > threshold(3))) = 3;
	InputSignal(find(InputSignal > threshold(2) & InputSignal <= threshold(3))) = 2;
	InputSignal(find(InputSignal > threshold(1) & InputSignal <= threshold(2))) = 1;
	InputSignal(find(InputSignal <= threshold(1))) = 0;
	
	%% Error Counting
	SymErrorNum = length(find(InputSignal ~= OriginalData));
	% The bit error number is 1/sym when only 1 bit changes and 2/sym when both bits change.
	BitErrorNum = SymErrorNum + length(find((InputSignal == 3) & (OriginalData == 0))) ...
														+ length(find((InputSignal == 2) & (OriginalData == 1))) ...
														+ length(find((InputSignal == 1) & (OriginalData == 2))) ...
														+ length(find((InputSignal == 0) & (OriginalData == 3)));
	SymErrorRate = SymErrorNum / length(OriginalData);
	BitErrorRate = BitErrorNum / (2 * length(OriginalData));
	