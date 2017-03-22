function [output, w, costs] = lmsFFEqualize(InputSignal, TrainingSignal, varargin)
	% This function performs the feed forward equalization with LMS algorithm.
	% For now, the training will use all the %InputSignal% and will be performed
	% %epoch% times. The %alpha% is the learning rate, which should be chosen
	% carefully with the help of curve of convergence. After training the equalization
	% will be performed and output the result.
	%
	% input: 
	%     InputSignal
	%       The input signal to be equalized.
	%     TrainingSignal
	%       The actual signal to be equalized to.
	%     FFETaps (optional)
	%       The numbers of FFE taps which must be odd.
	%       Default: 5
	%     alpha (optional)
	%       The learning rate of LMS algorithm.
	%       Default: 0.01
	%     epoch (optional)
	%       The epoch of the learning of LMS through all the input signal.
	%       Default: 1
	% output:
	%     output
	%       The equalized signal with the same length of InputSignal
	%       Size: length(InputSignal), 1
	%     w
	%       Weights of FFE
	%       Size: FFETaps, 1
	%     costs
	%       The costs after each training epoch, which is used to draw a 
	%       curve of convergence and thus determine the best learning rate.
	
	%% Parameter Checking
	narginchk(2, 5);
	
	if nargin == 2
		FFETaps = 5;
	else
		FFETaps = varargin{1};
	end
	if nargin <= 3
		alpha = 0.01;
	else
		alpha = varargin{2};
	end
	if nargin <= 4
		epoch = 1;
	else
		epoch = varargin{3};
	end
	
	% FFETaps must equals to a odd number
	if mod(FFETaps, 2) == 0
		error('lmsFeedForwardEqualize:argChk', 'FFE taps must be odd');
	end
	
	%% Signal Normalization and Duplication
	% InputSignal and TrainingSignal is normalized to the range between 0-1.
	InputSignal = InputSignal - min(InputSignal);
	InputSignal = InputSignal / max(InputSignal);
	TrainingSignal = TrainingSignal - min(TrainingSignal);
	TrainingSignal = TrainingSignal / max(TrainingSignal);
	
	% Both signal is duplicated for better performance
	InputSignalDup = repmat(InputSignal, 2, 1);
	TrainingSignalDup = repmat(TrainingSignal, 2, 1);
	% Zero Padding for input signal
	InputSignalZP = [zeros(floor(FFETaps/2), 1); InputSignalDup; zeros(floor(FFETaps/2), 1)];
	
	%% Weights Initializing
	w = zeros(FFETaps, 1);
	w(floor(length(w)/2) + 1) = 1;
	
	%% Training epoch times
	costs = zeros(epoch, 1);
	for n = 1 : epoch
		for i = 1 : length(InputSignalZP) - FFETaps + 1
			y(i) = w' * InputSignalZP(i : i + FFETaps - 1);
			w = w - alpha * (y(i) - TrainingSignalDup(i)) * InputSignalZP(i : i + FFETaps - 1);
			costs(n) = costs(n) + 0.5 * ((y(i) - TrainingSignalDup(i)) ^ 2);
		end
		% Record the cost/error of each epoch
		costs(n) = costs(n) / (length(InputSignalZP) - FFETaps + 1);
	end
	
	%% Using Trained Weights to Equalize Data
	for i = 1 : length(InputSignalZP) - FFETaps + 1
		y(i) = w' * InputSignalZP(i : i + FFETaps - 1);
	end
	
	y = y';
	
	% TODO choose a half of the output
	output = y(1 : length(y) / 2);
	% output = y(length(y) / 2 + 1 : end);
