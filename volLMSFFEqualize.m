function [output, w, costs] = volLMSFFEqualize(InputSignal, TrainingSignal, varargin)
	% This function performs the volterra feedforward equalization with LMS algorithm.
	% For now only 1st-3rd order are supported. 1st order must be included, while 2nd 
	% and 3rd orders are optional and controlled by the %en2ndOrder% and %en3rdOrder%
	% flags. The learning rate of different order can be different by adjusting the 
	% %alpha1st%, %alpha2nd% and %alpha3rd%. The equalizer will be trained on the 
	% %InputSignal% %epoch% times and will then perform a equalization.
	% 
	% input:
	%     InputSignal
	%       The input signal to be equalized.
	%     TrainingSignal
	%       The actual signal to be equalized to.
	%     chanLen (optional)
	%       The channel length. Must be odd.
	%       Default: 5
	%     alpha1st (optional)
	%       The learning rate of 1st-order kernel LMS algorithm.
	%       Default: 0.01
	%     epoch (optional)
	%       The epoch of the learning through all the input signal.
	%       Default: 1
	%     en2ndOrder (optional)
	%       The Flag of whether 2nd-order kernel is enabled.
	%       Default: true
	%     en3rdOrder (optional)
	%       The Flag of whether 3rd-order kernel is enabled.
	%       Default: true
	%     alpha2nd (optional)
	%       The learning rate of 2nd-order kernel LMS algorithm.
	%       Default: alpha1st
	%     alpha3rd (optional)
	%       The learning rate of 3rd-order kernel LMS algorithm.
	%       Default: alpha1st
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
	
	% TODO test
	
	%% Parameter Checking
	narginchk(2, 9);
	
	if nargin == 2
		chanLen = 5;
	else
		chanLen = varargin{1};
	end
	if nargin <= 3
		alpha1st = 0.01;
	else
		alpha1st = varargin{2};
	end
	if nargin <= 4
		epoch = 1;
	else
		epoch = varargin{3};
	end
	if nargin <= 5
		en2ndOrder = true;
	else
		en2ndOrder = varargin{4};
	end
	if nargin <= 6
		en3rdOrder = true;
	else
		en3rdOrder = varargin{5};
	end
	if nargin <= 7
		alpha2nd = alpha1st;
	else
		alpha2nd = varargin{6};
	end
	if nargin <= 8
		alpha3rd = alpha1st;
	else
		alpha3rd = varargin{7};
	end
	
	% channel length must equals to a odd number
	if mod(chanLen, 2) == 0
		error('volLmsFeedForwardEqualize:argChk', 'Channel length must be odd');
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
	
	%% Calculating the kernel/weight size
	weightSize = chanLen;
	if en2ndOrder == true
		kernel2ndSize = 0;
		for k = 1 : chanLen
			for m = k : chanLen
				kernel2ndSize = kernel2ndSize + 1;
			end
		end
		weightSize = weightSize + kernel2ndSize;
	end
	if en3rdOrder == true
		kernel3rdSize = 0;
		for k = 1 : chanLen
			for m = k : chanLen
				for n = m : chanLen
					kernel3rdSize = kernel3rdSize + 1;
				end
			end
		end
		weightSize = weightSize + kernel3rdSize;
	end
	% Define weights vector
	w = zeros(weightSize, 1);
	w(ceil(chanLen / 2)) = 1;
	
	%% Training epoch times
	costs = zeros(epoch, 1);
	for n = 1 : epoch
		for i = 1 : length(InputSignalZP) - chanLen + 1
			% Forming 2nd and 3rd kernel
			x = InputSignalZP(i : i + chanLen - 1);
			kernel = x;
			if en2ndOrder == true
				t = 0;
				for k = 1 : chanLen
					for m = k : chanLen
						t = t + 1;
						kernel2nd(t, :) = x(k) * x(m);
					end
				end
				kernel = [kernel; kernel2nd];
			end
			if en3rdOrder == true
				t = 0;
				for k = 1 : chanLen
					for m = k : chanLen
						for n = m : chanLen
							t = t + 1;
							kernel3rd(t, :) = x(k) * x(m) * x(n);
						end
					end
				end
				kernel = [kernel; kernel3rd];
			end
			% make a equalization
			y(i) = w' * kernel;
			% Construct a diagnose matrix for different learning rate of different order
			alpha = repmat(alpha1st, chanLen, 1);
			if en2ndOrder == true
				alpha = [alpha; repmat(alpha2nd, kernel2ndSize, 1)];
			end
			if en3rdOrder == true
				alpha = [alpha; repmat(alpha3rd, kernel3rdSize, 1)];
			end
			alpha = diag(alpha);
			% learning step
			w = w - alpha * (y(i) - TrainingSignalDup(i)) * kernel;
			costs(n) = costs(n) + 0.5 * ((y(i) - TrainingSignalDup(i)) ^ 2);
		end
		% Record the cost/error of each epoch
		costs(n) = costs(n) / (length(InputSignalZP) - chanLen + 1);
	end

	%% Using Trained Weights to Equalize Data
	for i = 1 : length(InputSignalZP) - chanLen + 1
		% Forming 2nd and 3rd kernel
		x = InputSignalZP(i : i + chanLen - 1);
		kernel = x;
		if en2ndOrder == true
			t = 0;
			for k = 1 : chanLen
				for m = k : chanLen
					t = t + 1;
					kernel2nd(t, :) = x(k) * x(m);
				end
			end
			kernel = [kernel; kernel2nd];
		end
		if en3rdOrder == true
			t = 0;
			for k = 1 : chanLen
				for m = k : chanLen
					for n = m : chanLen
						t = t + 1;
						kernel3rd(t, :) = x(k) * x(m) * x(n);
					end
				end
			end
			kernel = [kernel; kernel3rd];
		end
		% make a equalization
		y(i) = w' * kernel;
	end
	
	y = y';
	
	% TODO choose a half of the output
	output = y(1 : length(y) / 2);
	% output = y(length(y) / 2 + 1 : end);
	