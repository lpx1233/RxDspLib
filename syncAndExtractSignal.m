function [ExtractedSignalUS, OriginalSignalUS] = syncAndExtractSignal(SampledData, OriginalData, OverSamplingRatio, UpSamplingRatio)
	% This function performs the synchronization and signal extraction for the 
	% sampled signal from DSO. The %SampledData% will be down sampled by the 
	% %OverSamplingRatio%/%UpSamplingRatio% and correlated with the %OriginalData%
	% upsampled by %UpSamplingRatio%. Then the index of max correlating coeffients 
	% will be found and the synchronized data will be extracted.
	%
	% input:
	% 		SampledData
	%				The data sampled by DSO.
	%		  OriginalData
	%				The data which is transmitted with 1 symbol per second.
	%		  OverSamplingRatio
	%				The ratio of down sampling progress, which is usually the DSO sample 
	% 			rate dividing signal symbol rate
	%		  UpSamplingRatio (optional)
	%				The ratio of up sampling progress.
	%				Default value: 1
	% output:
	%     ExtractedSignalUS
	%       The synchronized and extracted signal upsampled by %UpSamplingRatio%.
	%       Size: length(OriginalData)*UpSamplingRatio, 1
	%     OriginalSignalUS
	%       The original signal upsampled by %UpSamplingRatio%.
	%       Size: length(OriginalData)*UpSamplingRatio, 1
	
	% Parameters checking
	narginchk(3,4);

	if ~exist('UpSamplingRatio','var') || isempty(UpSamplingRatio)
		UpSamplingRatio = 1;
	end
	
	% Downsampling
	DownSampledData = SampledData(1:OverSamplingRatio/UpSamplingRatio:end, 1);
	
	% Preparing transmitted data
	OriginalDataRemapped = (OriginalData - 1.5) * 2;
	OriginalSignalUS = reshape(repmat(OriginalDataRemapped, 1, UpSamplingRatio)', UpSamplingRatio * numel(OriginalDataRemapped), 1);
	
	% Correlation
	CorrelationResult = conv(DownSampledData(1:end), conj(OriginalSignalUS(end:-1:1)));
	[MaxCorr, index] = max(CorrelationResult);
	
	ExtractedSignalUS = DownSampledData(index-length(OriginalData)*UpSamplingRatio+1 : index);
