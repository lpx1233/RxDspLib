function [output, ChnlCoeffs, costs] = mlseEqualize(InputSignal, TrainingSignal, ...
																										TraceBackLen, SamplesPerSymbol, ...
																										ChanLen, ModFormat, ...
																										epoch, alpha)
	% This function performs equalization based on maximum likelhood sequence estimation.
	% First, the InputSignal and TrainingSignal will be normalized to 0-1.
	% Then, the inverse linear FFE, which is a channel estimation, is performed
	% using a ChanLen taps FFE with a learning rate of alpha and training epoch times.
	% After training, the MLSE equalization using MATLAB Comm. Toolbox will be performed 
	% and then the result will be returned. The parameter of MLSE equalizer such as TraceBackLen
	% and SamplesPerSymbol can be specified. The ModFormat determine the signal constellation map.
	%
	% input: 
	%     InputSignal
	%       The input signal to be equalized.
	%     TrainingSignal
	%       The actual signal to be equalized to.
	%     TraceBackLen (optional)
	%       The trace back length of mlse equalizer
	%       Default: 4
	%     SamplesPerSymbol (optional)
	%       The numbers of samples per symbol.
	%       Default: 1
	%     ChanLen (optional)
	%       The Length of channel coefficients.
	%       Default: 5
	%     ModFormat (optional)
	%       The modulation format: 'PAM4' or 'NRZ'. More modulation format will be added.
	%       Default: 'PAM4'
	%     epoch (optional)
	%       The number of epochs to running in the channel estimation step.
	%       Default: 5
	%     alpha (optional)
	%       The learning rate of LMS algorithm for channel estimation.
	%       Default: 0.01
	% output:
	%     output
	%       The equalized signal with the same length of InputSignal
	%       Size: length(InputSignal), 1
	%     ChnlCoeffs
	%       Channel coefficients
	%       Size: ChanLen, 1
	%     costs
	%       The costs of channel estimation after each training epoch, which is 
	%       used to draw a curve of convergence and thus determine the best learning rate.
	
	%% Parameter Checking
	narginchk(2, 8);
	
	if ~exist('TraceBackLen','var') || isempty(TraceBackLen)
		TraceBackLen = 4;
	end
	if ~exist('SamplesPerSymbol','var') || isempty(SamplesPerSymbol)
		SamplesPerSymbol = 1;
	end
	if ~exist('ChanLen','var') || isempty(ChanLen)
		ChanLen = 5;
	end
	if ~exist('ModFormat','var') || isempty(ModFormat)
		ModFormat = 'PAM4';
	end
	if ~exist('epoch','var') || isempty(epoch)
		epoch = 5;
	end
	if ~exist('alpha','var') || isempty(alpha)
		alpha = 0.01;
	end
	
	%% Signal Normalization and Duplication
	% InputSignal and TrainingSignal is normalized to the range between 0-1.
	InputSignal = InputSignal - min(InputSignal);
	InputSignal = InputSignal / max(InputSignal);
	TrainingSignal = TrainingSignal - min(TrainingSignal);
	TrainingSignal = TrainingSignal / max(TrainingSignal);
	
	%% Channel Estimation
	% Zero Padding for training signal
	TrainingSignalZP = [zeros(ChanLen - 1, 1); TrainingSignal];
	
	% Weights Initializing
	w = zeros(ChanLen, 1);
	w(ChanLen) = 1;
	
	% Training using LMS learning algorithm
	costs = zeros(epoch, 1);
	for n = 1 : epoch
		for i = 1 : length(TrainingSignalZP) - ChanLen + 1
			y(i) = w' * TrainingSignalZP(i : i + ChanLen - 1);
			w = w - alpha * (y(i) - InputSignal(i)) * TrainingSignalZP(i : i + ChanLen - 1);
			costs(n) = costs(n) + 0.5 * ((y(i) - InputSignal(i)) ^ 2);
		end
		% Record the cost/error of each epoch
		costs(n) = costs(n) / (length(TrainingSignalZP) - ChanLen + 1);
	end
	
	% Reverse w to get channel coefficients
	ChnlCoeffs = w(end:-1:1);

	%% MLSE equalization using MATLAB Commmunication Toolbox
	% Generating constellation map for different modulation format
	if ModFormat == 'PAM4'
		% const = [1/8; 3/8; 5/8; 7/8];
		const = [0; 1/3; 2/3; 1];
	elseif ModFormat == 'NRZ'
		const = [1/4; 3/4];
	end
	% Performing the mlse equalization
	output = mlseeq(InputSignal, ChnlCoeffs, const, TraceBackLen, 'rst', SamplesPerSymbol);
