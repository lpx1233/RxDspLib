clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));

DataLength = 500000;
NumOfSeq = 5;
rng('shuffle');
PAM4Data = randi(4, DataLength, 1) - 1;
PAM4Data = reshape(PAM4Data, DataLength / NumOfSeq, NumOfSeq);
PAM8Data = randi(8, DataLength, 1) - 1;
PAM8Data = reshape(PAM8Data, DataLength / NumOfSeq, NumOfSeq);
PAM16Data = randi(16, DataLength, 1) - 1;
PAM16Data = reshape(PAM16Data, DataLength / NumOfSeq, NumOfSeq);

if ~exist('Original Data', 'dir')
    mkdir('Original Data');
end

for i = 1 : NumOfSeq
    for j = [4, 8, 16]
        csvwrite(['Original Data\pam', num2str(j), '_', num2str(i), '.csv'], eval(['PAM', num2str(j), 'Data(:, ', num2str(i), ')']));
    end
end

PAM4Data = PAM4Data / 3;
PAM8Data = PAM8Data / 7;
PAM16Data = PAM16Data / 15;
