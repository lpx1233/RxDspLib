function OriginalData = generateData(OriginalDataLength, PAM4Flag, ...
																		 NewPRBSGenerationFlag, SyncZerosLength)
	% This function generate NRZ or PAM4 data from PRBS generator or local file 
	% and then generate the data with synchronization header which is then 
	% loaded to PPG.
	% 
	% input:(input parameter are all optional, the order of args follows the order below)
	%		  OriginalDataLength
	%				The length of data to be generated.
	%				Default value: 2^12
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
	%		  SyncZerosLength
	%				The length of zeros in the synchronization header for data to PPG.
	%				Default value: 50
	% output: 
	%     OriginalData
	%       The original PAM4 or NRZ data without synchronization header in row vector.
	%       Size: OriginalDataLength, 1
	
	narginchk(0,4);
	
	if ~exist('OriginalDataLength','var') || isempty(OriginalDataLength)
		OriginalDataLength = 2^12;
	end
	
	if ~exist('PAM4Flag','var') || isempty(PAM4Flag)
		PAM4Flag = 1;
	end
	
	if ~exist('NewPRBSGenerationFlag','var') || isempty(NewPRBSGenerationFlag)
		NewPRBSGenerationFlag = 0;
	end
	
	if ~exist('SyncZerosLength','var') || isempty(SyncZerosLength)
		SyncZerosLength = 50;
	end
	
	% change the current directory to the folder which contains this m file
	cd(fileparts(which(mfilename)));
	
	%% define parameter
	PathToOriginalData = '.\Original Data\'; % relative path
	OriginalDataFileName = strcat('Original_Data_', num2str(OriginalDataLength));
	
	if exist('.\Original Data', 'dir') == 0
		mkdir('Original Data');
	end
	
	%% PRBS generation
	if (NewPRBSGenerationFlag == 1) || (exist(strcat(PathToOriginalData, OriginalDataFileName, '.txt'), 'file') == 0)
		h = commsrc.pattern('SamplingFrequency', 10000, ...
							'SamplesPerSymbol', 1, ...
							'PulseType', 'NRZ', ...
							'OutputLevels', [0 1], ...
							'RiseTime', 0, ...
							'FallTime', 0, ...
							'DataPattern', 'PRBS12'); % set the prbs pattern
		OriginalData = generate(h, OriginalDataLength);
		OriginalData = double(xor(1, OriginalData));
		
		% save the prbs data to .\Original Data\Original_Data_4096.txt
		fid = fopen(strcat(PathToOriginalData, OriginalDataFileName, '.txt'), 'w'); 
		fprintf(fid, '%d\r\n', OriginalData); 
		fclose(fid);
		
		% save the data with sync header to PPG to .\Original Data\Data2PPG.txt
		fid = fopen(strcat(PathToOriginalData, 'Data2PPG', '.txt'), 'w');
		fprintf(fid, '%d\r\n', [zeros(SyncZerosLength, 1); OriginalData]);
		fclose(fid);
	end
	OriginalData = importdata(strcat(PathToOriginalData, OriginalDataFileName, '.txt'));
	
	%% PAM4 generation
	% use 2 same prbs data, delay one with 2 sym, multiply the other one with 2
	% and then add the 2 data to form PAM4 data
	if PAM4Flag == 1
		OriginalData_port1 = OriginalData;
		OriginalData_port2 = [1; 1; ~(OriginalData(1 : length(OriginalData)-2))];
		OriginalData = 2*OriginalData_port1 + OriginalData_port2;
		% save the PAM4 data to .\Original Data\Original_Data_4096_PAM4.txt
		fid = fopen(strcat(PathToOriginalData, OriginalDataFileName, '_PAM4.txt'), 'w');
		fprintf(fid, '%d\r\n', OriginalData);
		fclose(fid);
	end
