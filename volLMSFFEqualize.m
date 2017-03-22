function [output, w, costs] = volLMSFFEqualize(InputSignal, TrainingSignal, varargin)
	% This function performs the volterra feedforward equalization with LMS algorithm.
	% 
	% input:
	%     InputSignal
	%       The input signal to be equalized.
	%     TrainingSignal
	%       The actual signal to be equalized to.
	%     chanLen (optional)
	%       The channel length. Must be odd.
	%       Default: 5
	%     epoch
	%       The epoch of the learning through all the input signal.
	%       Default: 1
	%
	% output:
	%     output
	%       The equalized signal with the same length of InputSignal
	%       Size: length(InputSignal), 1
	%     w
	%       Weights of volterra series
	%       Size: chanLen + kernel2ndSize + kernel3rdSize, 1
	%     costs
	%       The costs after each training epoch, which is used to draw a 
	%       curve of convergence and thus determine the best learning rate.
	
	%% Parameter Checking
	narginchk(2, 4);
	
	if nargin == 2
		chanLen = 5;
	else
		chanLen = varargin{1};
	end
	if nargin <= 3
		epoch = 1;
	else
		epoch = varargin{2};
	end
	
	% FFETaps must equals to a odd number
	if mod(chanLen, 2) == 0
		error('lmsFeedForwardEqualize:argChk', 'Channel length must be odd');
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
	InputSignalZP = [zeros(floor(chanLen/2), 1); InputSignalDup; zeros(floor(chanLen/2), 1)];
	
	%% Calculating the kernel size
	kernel2ndSize = 0;
	for k = 1 : chanLen
		for m = k : chanLen
			kernel2ndSize = kernel2ndSize + 1;
		end
	end
	kernel3rdSize = 0;
	for k = 1 : chanLen
		for m = k : chanLen
			for n = m : chanLen
				kernel3rdSize = kernel3rdSize + 1;
			end
		end
	end
	
	% Define weights vector
	w = zeros(chanLen + kernel2ndSize + kernel3rdSize, 1);
	w(ceil(chanLen / 2)) = 1;
	
	%% Training epoch times
	costs = zeros(epoch, 1);
	for n = 1 : epoch
		for i = 1 : length(InputSignalZP) - chanLen + 1
			%% Forming 2nd and 3rd kernel
			x = InputSignalZP(i : i + chanLen - 1);
			t = 0;
			for k = 1 : chanLen
				for m = k : chanLen
					t = t + 1;
					kernel2nd(t, :) = x(k) * x(m);
				end
			end
			t = 0;
			for k = 1 : chanLen
				for m = k : chanLen
					for n = m : chanLen
						t = t + 1;
						kernel3rd(t, :) = x(k) * x(m) * x(n);
					end
				end
			end
			kernel = [InputSignalZP(i : i + chanLen - 1); kernel2nd; kernel3rd];
			y(i) = w' * kernel;
			
			% TODO modify the w-updating algorithm and costs
			w = w - alpha * (y(i) - TrainingSignalDup(i)) * InputSignalZP(i : i + FFETaps - 1);
			costs(n) = costs(n) + 0.5 * ((y(i) - TrainingSignalDup(i)) ^ 2);
		end
		% Record the cost/error of each epoch
		costs(n) = costs(n) / (length(InputSignalZP) - FFETaps + 1);
	end

	%% Using Trained Weights to Equalize Data
	for i = 1 : length(InputSignalZP) - chanLen + 1
		%% Forming 2nd and 3rd kernel
		x = InputSignalZP(i : i + chanLen - 1);
		t = 0;
		for k = 1 : chanLen
			for m = k : chanLen
				t = t + 1;
				kernel2nd(t, :) = x(k) * x(m);
			end
		end
		t = 0;
		for k = 1 : chanLen
			for m = k : chanLen
				for n = m : chanLen
					t = t + 1;
					kernel3rd(t, :) = x(k) * x(m) * x(n);
				end
			end
		end
		kernel = [InputSignalZP(i : i + chanLen - 1); kernel2nd; kernel3rd];
		y(i) = w' * kernel;
	end
	