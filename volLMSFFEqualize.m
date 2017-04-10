function [output, w, costs] = volLMSFFEqualize(InputSignal, TrainingSignal, ...
																								AlgType, epoch, ...
																								ChanLen1st, Alpha1st, ...
																								ChanLen2nd, Alpha2nd, ...
																								ChanLen3rd, Alpha3rd)
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
	%     AlgType
	%       'lms' for LMS or 'rls' for RLS.
	%     epoch (optional)
	%       The epoch of the learning through all the input signal.
	%       Default: 1
	%     ChanLen1st (optional)
	%       The channel length of 1st order kernel.
	%       Default: 5
	%     Alpha1st (optional)
	%       The 1st order kernel learning rate of LMS algorithm or the forgetting factor of RLS.
	%       Default: 0.01 for AlgType = 'lms', 0.99 for AlgType = 'rls'
	%     ChanLen2nd (optional)
	%       The channel length of 2nd order kernel.
	%       Default: ChanLen1st
	%     Alpha2nd (optional)
	%       The 2nd order kernel learning rate of LMS algorithm or the forgetting factor of RLS.
	%       In rls algorithm, this parameter is forced to be Alpha1st.
	%       Default: Alpha1st
	%     ChanLen3rd (optional)
	%       The channel length of 3rd order kernel.
	%       Default: ChanLen1st
	%     Alpha3rd (optional)
	%       The 3rd order kernel learning rate of LMS algorithm or the forgetting factor of RLS.
	%       In rls algorithm, this parameter is forced to be Alpha1st.
	%       Default: Alpha1st
	%
	% output:
	%     output
	%       The equalized signal with the same length of InputSignal
	%       Size: length(InputSignal), 1
	%     w
	%       Weights of volterra series
	%       Size: chanLen + Kernel2ndSize + Kernel3rdSize, 1
	%     costs
	%       The costs after each training epoch, which is used to draw a 
	%       curve of convergence and thus dete_raermine the best learning rate.
	
	% TODO test
	
	%% Paramete_raer Checking
	narginchk(2, 10);
	
	% algorithm type must be lms or rls
	if (AlgType ~= 'lms') & (AlgType ~= 'rls')
		error('volterraFeedForwardEqualize:argChk', 'AlgType must be lms or rls');
	end
	
	if ~exist('epoch','var') || isempty(epoch)
		epoch = 1;
	end
	
	if ~exist('ChanLen1st','var') || isempty(ChanLen1st)
		ChanLen1st = 5;
	end
	if ~exist('Alpha1st','var') || isempty(Alpha1st)
		if AlgType == 'lms'
			Alpha1st = 0.01;
		elseif AlgType == 'rls'
			Alpha1st = 0.99;
		else
			error('volterraFeedForwardEqualize:argChk', 'AlgType must be lms or rls');
		end
	end
	
	if ~exist('ChanLen2nd','var') || isempty(ChanLen2nd)
		ChanLen2nd = ChanLen1st;
	end
	if ~exist('Alpha2nd','var') || isempty(Alpha2nd)
		Alpha2nd = Alpha1st;
	end
	
	if ~exist('ChanLen3rd','var') || isempty(ChanLen3rd)
		ChanLen3rd = ChanLen1st;
	end
	if ~exist('Alpha3rd','var') || isempty(Alpha3rd)
		ChanLen3rd = Alpha1st;
	end

	% channel length must equals to a odd number
	if (mod(ChanLen1st, 2) == 0) | (mod(ChanLen2nd, 2) == 0) | (mod(ChanLen3rd, 2) == 0)
		error('volterraFeedForwardEqualize:argChk', 'Channel length must be odd');
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
	MaxChanLen = max([ChanLen1st, ChanLen2nd, ChanLen3rd]);
	InputSignalZP = [zeros(floor(MaxChanLen/2), 1); InputSignalDup; zeros(floor(MaxChanLen/2), 1)];
	
	%% Calculating the kernel size
	KernelSize = ChanLen1st;
	if ChanLen2nd ~= 0
		Kernel2ndSize = 0;
		for k = 1 : ChanLen2nd
			for m = k : ChanLen2nd
				Kernel2ndSize = Kernel2ndSize + 1;
			end
		end
		KernelSize = KernelSize + Kernel2ndSize;
	end
	if ChanLen3rd ~= 0
		Kernel3rdSize = 0;
		for k = 1 : ChanLen3rd
			for m = k : ChanLen3rd
				for n = m : ChanLen3rd
					Kernel3rdSize = Kernel3rdSize + 1;
				end
			end
		end
		KernelSize = KernelSize + Kernel3rdSize;
	end
	
	%% Weights Searching using RAE
	% Generating input matrix
	X = zeros(length(InputSignalDup), KernelSize);
	for i = 1 : ChanLen1st
		X(:, i) = InputSignalZP(floor(MaxChanLen/2) + i - floor(ChanLen1st/2) : ...
														floor(MaxChanLen/2) + i - floor(ChanLen1st/2) + length(InputSignalDup) - 1);
	end
	t = ChanLen1st;
	if ChanLen2nd ~= 0
		for k = 1 : ChanLen2nd
			for m = k : ChanLen2nd
				t = t + 1;
				X(:, t) = InputSignalZP(floor(MaxChanLen/2) + k - floor(ChanLen2nd/2) : ...
														floor(MaxChanLen/2) + k - floor(ChanLen2nd/2) + length(InputSignalDup) - 1) .*...
									InputSignalZP(floor(MaxChanLen/2) + m - floor(ChanLen2nd/2) : ...
														floor(MaxChanLen/2) + m - floor(ChanLen2nd/2) + length(InputSignalDup) - 1);
			end
		end
	end
	if ChanLen3rd ~= 0
		for k = 1 : ChanLen3rd
			for m = k : ChanLen3rd
				for n = m : ChanLen3rd
					t = t + 1;
					X(:, t) = InputSignalZP(floor(MaxChanLen/2) + k - floor(ChanLen3rd/2) : ...
															floor(MaxChanLen/2) + k - floor(ChanLen3rd/2) + length(InputSignalDup) - 1) .*...
										InputSignalZP(floor(MaxChanLen/2) + m - floor(ChanLen3rd/2) : ...
															floor(MaxChanLen/2) + m - floor(ChanLen3rd/2) + length(InputSignalDup) - 1) .*...
										InputSignalZP(floor(MaxChanLen/2) + n - floor(ChanLen3rd/2) : ...
															floor(MaxChanLen/2) + n - floor(ChanLen3rd/2) + length(InputSignalDup) - 1);
				end
			end
		end
	end
	
	% RAE process: Q is the output selected kernel matrix
	y = TrainingSignalDup;
	p = 1;
	epsilon = 0.01;
	ete_rae = zeros(KernelSize, KernelSize);
	
	for i = 1 : KernelSize
		w_rae = (X(:, i)' * y) / (X(:, i)' * X(:, i));
		ete_rae(1, i) = y' * y - w_rae * X(:, i)' * y;
	end
	[ete_rae_min, i] = min(ete_rae(1, :));
	w_rae = (X(:, i)' * y) / (X(:, i)' * X(:, i));
	Q = [X(:, i)];
	if i <= ChanLen1st
		KernelOrder = [1];
	elseif i <= ChanLen1st + Kernel2ndSize
		KernelOrder = [2];
	else
		KernelOrder = [3];
	end
	R = Q' / (Q' * Q);
	p = p + 1;

	while (ete_rae_min > (epsilon * length(InputSignalDup))) | (p ~= KernelSize)
		for i = 1 : KernelSize
			if ismember(X(:, i)', Q', 'rows') == 1
				continue
			end
			z_rae = R * X(:, i);
			delta_rae = 1 / (X(:, i)' * X(:, i) - X(:, i)' * Q * z_rae);
			g_rae = delta_rae * (z_rae' * Q' - X(:, i)');
			w_rae_temp = [w_rae + z_rae * g_rae * y; -g_rae * y];
			ete_rae(p, i) = y' * y - y' * [Q, X(:, i)] * w_rae_temp;
		end
		[ete_rae_min, i] = min(ete_rae(p, :));
		z_rae = R * X(:, i);
		delta_rae = 1 / (X(:, i)' * X(:, i) - X(:, i)' * Q * z_rae);
		g_rae = delta_rae * (z_rae' * Q' - X(:, i)');
		w_rae = [w_rae + z_rae * g_rae * y; -g_rae * y];
		Q = [Q, X(:, i)];
		if i <= ChanLen1st
			KernelOrder = [KernelOrder; 1];
		elseif i <= ChanLen1st + Kernel2ndSize
			KernelOrder = [KernelOrder; 2];
		else
			KernelOrder = [KernelOrder; 3];
		end
		R = [R + z_rae * g_rae; -g_rae];
		p = p + 1;
	end
	
	% Define and randomly init weights vector to (-1, 1)
	w = rand(length(Q, 2), 1);
	w = 2 * w - 1;
	
	%% Training epoch times
	costs = zeros(epoch, 1);
	y = zeros(size(TrainingSignalDup));
	if AlgType == 'lms'
		for iter = 1 : epoch
			for i = 1 : length(Q, 1)
				y(i) = Q(i, :) * w;
				% Construct a diagnose matrix for different learning rate of different order
				alpha = [alpha1st; alpha2nd; alpha3rd];
				alpha = diag(alpha(KernelOrder));
				% learning step
				w = w - alpha * (y(i) - TrainingSignalDup(i)) * kernel;
				costs(iter) = costs(iter) + 0.5 * ((y(i) - TrainingSignalDup(i)) ^ 2);
			end
			% Record the cost/error of each epoch
			costs(iter) = costs(iter) / (length(InputSignalZP) - chanLen + 1);
		end
	elseif AlgType == 'rls'
		% TODO implement rls algorithm
		% The RLS learning algorithm
		Sd = eye(length(Q, 2));
		for n = 1 : epoch
			for i = 1 : length(InputSignalZP) - FFETaps + 1
				x = InputSignalZP(i : i + FFETaps - 1);
				e = TrainingSignalDup(i) - w' * x;
				phi = Sd * x;
				Sd = (1 / Alpha1st) * (Sd - (phi * phi') / (Alpha1st + phi' * x));
				w = w + e * Sd * x;
				costs(n) = costs(n) + 0.5 * (e ^ 2);
			end
			% Record the cost/error of each epoch
			costs(n) = costs(n) / (length(InputSignalZP) - FFETaps + 1);
		end
	end
	
	%% Using Trained Weights to Equalize Data
	% TODO modify these code to be functional
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
	