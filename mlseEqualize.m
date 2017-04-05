function [output, ChnlCoeffs, costs] = mlseEqualize(InputSignal, TrainingSignal, varargin)
	% This function performs the feed forward equalization with LMS or RLS algorithm.
	% First, the %InputSignal% and %TrainingSignal% will be normalized to 0-1.
	% The training will use all the %InputSignal% and will be performed %epoch% times.
	% The %alpha% is the learning rate of LMS or the forgetting factor of RLS, 
	% which should be chosen carefully with the help of curve of convergence. 
	% After training, the equalization will be performed and then the result will 
	% be returned.
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
	
	if nargin == 2
		TraceBackLen = 4;
	else
		TraceBackLen = varargin{1};
	end
	if nargin <= 3
		SamplesPerSymbol = 1;
	else
		SamplesPerSymbol = varargin{2};
	end
	if nargin <= 4
		ChanLen = 5;
	else
		ChanLen = varargin{3};
	end
	if nargin <= 5
		ModFormat = 'PAM4';
	else
		ModFormat = varargin{4};
	end
	if nargin <= 6
		epoch = 5;
	else
		epoch = varargin{5};
	end
	if nargin <= 7
		alpha = 0.01;
	else
		alpha = varargin{6};
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
	% Generating constallation map for different modulation format
	if ModFormat == 'PAM4'
		% const = [1/8; 3/8; 5/8; 7/8];
		const = [0; 1/3; 2/3; 1];
	elseif ModFormat == 'NRZ'
		const = [1/4; 3/4];
	end
	% Performing the mlse equalization
	output = mlseeq(InputSignal, ChnlCoeffs, const, TraceBackLen, 'rst', SamplesPerSymbol);
