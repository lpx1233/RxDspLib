clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

if ~exist('Original Data', 'dir')
    error('Original Data fold does not exist');
end

DataLength = 500000;
NumOfSeq = 5;
PAM4Data = zeros(DataLength / NumOfSeq, NumOfSeq);
PAM8Data = zeros(DataLength / NumOfSeq, NumOfSeq);
PAM16Data = zeros(DataLength / NumOfSeq, NumOfSeq);

for i = 1 : NumOfSeq
  for j = [4, 8, 16]
    switch j
      case 4
        PAM4Data(:, i) = importdata(['.\Original Data\pam4_', num2str(i), '.csv']);
      case 8
        PAM8Data(:, i) = importdata(['.\Original Data\pam8_', num2str(i), '.csv']);
      case 16
        PAM16Data(:, i) = importdata(['.\Original Data\pam16_', num2str(i), '.csv']);
    end
  end
end

PAM4Data = PAM4Data / 3;
PAM8Data = PAM8Data / 7;
PAM16Data = PAM16Data / 15;
