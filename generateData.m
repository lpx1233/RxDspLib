function OriginalData = generateData(PAM4Flag, NewPRBSGenerationFlag)
	% This function generate NRZ or PAM4 data from PRBS generator or local file
	% and then generate the data with synchronization header which is then
	% loaded to PPG.
	%
	% input:(input parameter are all optional, the order of args follows the order below)
	%		  PAM4Flag
	%				The flag for PAM4 generation:
	%			 		1: enable PAM4 generation
	%			 		0: disable PAM4 generation and output NRZ signal
	%				Default value: 1
	%		  NewPRBSGenerationFlag
	%				The flag for regenerating prbs data:
	%			  	1: regenerating new data
	%			  	0: load prbs data from ./Original Data/Original_Data_4096.txt
	%				Note that if the ./Original Data/Original_Data_4096 doesn't exist,
	%				the program will automatically generate new data whether this flag
	%				is set to 1 or 0.
	%				Default value: 0
	% output:
	%     OriginalData
	%       The original PAM4 or NRZ data without synchronization header in row vector.
	%       Size: OriginalDataLength, 1

	narginchk(0,2);

	if ~exist('PAM4Flag','var') || isempty(PAM4Flag)
		PAM4Flag = 1;
	end

	if ~exist('NewPRBSGenerationFlag','var') || isempty(NewPRBSGenerationFlag)
		NewPRBSGenerationFlag = 0;
	end

	% change the current directory to the folder which contains this m file
	cd(fileparts(which(mfilename)));

	%% define parameter
	PathToOriginalData = '.\Original Data\'; % relative path
	OriginalDataFileName = 'Original_Data';

	if exist('.\Original Data', 'dir') == 0
		mkdir('Original Data');
	end

	%% PRBS generation
	if (NewPRBSGenerationFlag == 1) || (exist(strcat(PathToOriginalData, OriginalDataFileName, '.txt'), 'file') == 0)
		% h = commsrc.pattern('SamplingFrequency', 10000, ...
							% 'SamplesPerSymbol', 1, ...
							% 'PulseType', 'NRZ', ...
							% 'OutputLevels', [0 1], ...
							% 'RiseTime', 0, ...
							% 'FallTime', 0, ...
							% 'DataPattern', 'PRBS12'); % set the prbs pattern
		% OriginalData = generate(h, OriginalDataLength);
		% OriginalData = double(xor(1, OriginalData));
		% Generating PRBS15
		fbconnection = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1];
		n = length(fbconnection);
		N = 2 ^ n - 1;
		OriginalData = zeros(N, 1);
		register = [zeros(1, n - 1), 1];
		OriginalData(1) = register(n);
		for i = 2 : N
			newregister(1) = mod(sum(fbconnection .* register), 2);
			for j = 2 : n
					newregister(j) = register(j - 1);
			end
			register = newregister;
			OriginalData(i) = register(n);
		end

		% save the prbs data to .\Original Data\Original_Data.txt
		fid = fopen(strcat(PathToOriginalData, OriginalDataFileName, '.txt'), 'w');
		fprintf(fid, '%d\r\n', OriginalData);
		fclose(fid);
	end
	OriginalData = importdata(strcat(PathToOriginalData, OriginalDataFileName, '.txt'));

	%% PAM4 generation
	% use 2 same prbs data, delay one with 2 sym, multiply the other one with 2
	% and then add the 2 data to form PAM4 data
	if PAM4Flag == 1
		OriginalData_port1 = OriginalData;
		% Day 2 Setup
		% shiftnum = 6;
		% OriginalData_port2 = [~(OriginalData(end - shiftnum + 1 : end));
		% 											~(OriginalData(1 : end - shiftnum))];
		% OriginalData = OriginalData_port1 + 2 * OriginalData_port2;

		% 20170428 Setup
		shiftnum = 6;
		OriginalData_port2 = [~(OriginalData(shiftnum + 1 : end));
													~(OriginalData(1 : shiftnum))];
		OriginalData = 2 * OriginalData_port1 + OriginalData_port2;

		% save the PAM4 data to .\Original Data\Original_Data_4096_PAM4.txt
		fid = fopen(strcat(PathToOriginalData, OriginalDataFileName, '_PAM4.txt'), 'w');
		fprintf(fid, '%d\r\n', OriginalData);
		fclose(fid);
	end
