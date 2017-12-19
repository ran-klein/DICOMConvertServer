% DCSFinalizeData(DCSCurrentData) - finalizes conversion of a parsed dicom
% series and saves a .mat file
%
% See also: dicomConvertServer

% By Ran Klein August-2007
% University of Ottawa Heart Institute

% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************

function [vol, hdr, infodcm] = DCSFinalizeData(DCSCurrentData)

hdr = DCSCurrentData.hdr;
chklist = DCSCurrentData.chklist;
infodcm = DCSCurrentData.infodcm;

hdr.DownSampleFactor_xy_z = [1 1];
hdr.nframes = length(chklist.frame.vals);
hdr.nplanes = length(chklist.slice.vals);
if strcmpi(chklist.slice.mode,'Explicit Location')
	hdr.pix_mm_z = abs(nanmean(diff(chklist.slice.vals)));
end

% Ensure that dynamic arrays have nframe elements
if length(hdr.frame_start)~=hdr.nframes
	hdr.frame_start = nan(hdr.nframes,1);
end
if length(hdr.frame_len)~=hdr.nframes
	hdr.frame_len = nan(hdr.nframes,1);
end
if length(hdr.quant_dynamic)~=hdr.nframes
	hdr.quant_dynamic = nan(hdr.nframes,1);
end
if length(hdr.PETTotalCounts)~=hdr.nframes
	hdr.PETTotalCounts = nan(hdr.nframes,1);
end
if length(hdr.PrimaryPromptsCountsAccumulated)~=hdr.nframes
	hdr.PrimaryPromptsCountsAccumulated = nan(hdr.nframes,1);
end
if length(hdr.ScatterFractionFactor)~=hdr.nframes
	hdr.ScatterFractionFactor = nan(hdr.nframes,1);
end
if length(hdr.DeadTimeFactor)~=hdr.nframes
	hdr.DeadTimeFactor = nan(hdr.nframes,1);
end

vol1 = [];
prevfilei = 0;
for fi = 1:hdr.nframes
	for si = 1:hdr.nplanes
		if strcmpi(get(DCSlogmsg,'UserData'),'HaltRequest')
			set(DCSlogmsg,'UserData','Halted');
			return;
		end
		
		% which file is needed
		filei = chklist.table(si,fi,1);
		
		if filei~=prevfilei
			info = dicominfo(chklist.files{filei});
			% 	info = mydicominfo(chklist.files{filei},DCSVars.dict,DCSVars.codes);
			% read file image
			temp = squeeze(dicomread(chklist.files{filei}));
			if isempty(temp)
				error('EmptyImg:SaveTarget:DicomConertServer','The image in first file was empty')
			end
			prevfilei = filei;
		end
		
		% preacllocate image variable according to dimensions
		if fi==1 && si==1
			DCSlogmsg('Allocating memory for image');
			if isempty(hdr.xdim), hdr.xdim = size(temp,1); end
			if isempty(hdr.ydim), hdr.ydim = size(temp,1); end
			if ~rem(hdr.transverseRotation,180)
				t = hdr.xdim;		hdr.xdim = hdr.ydim;		hdr.ydim = t; % swap dimensions
			end
			while isempty(vol1)
				try
					if hdr.nframes>1
						vol = zeros(ceil(hdr.xdim/hdr.DownSampleFactor_xy_z(1)) ,ceil(hdr.ydim/hdr.DownSampleFactor_xy_z(1)) , ceil(hdr.nplanes/hdr.DownSampleFactor_xy_z(2)), hdr.nframes, 'int16');
					end
					vol1 = zeros(ceil(hdr.xdim/hdr.DownSampleFactor_xy_z(1)) ,ceil(hdr.ydim/hdr.DownSampleFactor_xy_z(1)) , ceil(hdr.nplanes/hdr.DownSampleFactor_xy_z(2)));
				catch
					% Maintain a limited aspect ratio between slice thickness
					% and pixel size
					aspectr = (hdr.pix_mm_z*hdr.DownSampleFactor_xy_z(2)) / (hdr.pix_mm_xy*hdr.DownSampleFactor_xy_z(1));
					if aspectr<1/2
						hdr.DownSampleFactor_xy_z(2) = hdr.DownSampleFactor_xy_z(2)*2;
					elseif aspectr > 4
						hdr.DownSampleFactor_xy_z(1) = hdr.DownSampleFactor_xy_z(1)*2;
					else
						hdr.DownSampleFactor_xy_z = hdr.DownSampleFactor_xy_z*2;
					end
				end
			end
			if any(hdr.DownSampleFactor_xy_z>1)
				DCSlogmsg(['  - * Warning! Due to memory limitation the image was down sampled (x and y factor = ' num2str(hdr.DownSampleFactor_xy_z(1)) ' and z factor = ' num2str(hdr.DownSampleFactor_xy_z(2)) ') *']);
				DCSlogmsg(['  - Original dimensions: ' num2str(hdr.xdim) ' ' num2str(hdr.ydim) ' ' num2str(hdr.nplanes) ' ' num2str(hdr.nframes)]);
			end
			DCSlogmsg('Adding Frames:');
		end
		
		DCSlogmsg(['Adding Frames: |' '*'*ones(1,round(20*(hdr.nplanes*(fi-1)+si)/hdr.nframes/hdr.nplanes)) '_'*ones(1,20-round(20*(hdr.nplanes*(fi-1)+si)/hdr.nframes/hdr.nplanes)) '|'],2);
		
		if ~isfield(info,'RescaleSlope')
			if isfield(info,'WindowWidth')
				info.RescaleSlope = info.WindowWidth/double(info.LargestImagePixelValue);
			elseif isfield(info,'Manufacturer') && strcmpi(info.Manufacturer,'GE MEDICAL SYSTEMS') && isfield(info,'Private_0011_103b')
				info.RescaleSlope = info.Private_0011_103b;
			else
				info.RescaleSlope = 1;
			end
		end
		if ~isfield(info,'RescaleIntercept')
			if isfield(info,'WindowCenter')
				info.RescaleIntercept = (info.WindowCenter-info.WindowWidth/2);
			else
				info.RescaleIntercept = 0;
			end
		end
		
		try
			ti=chklist.table(si,fi,2);
			if hdr.DownSampleFactor_xy_z(2) == 1
				vol1(:,:,si) = info.RescaleSlope*imrotate(imresize(double(temp(:,:,ti)),1/hdr.DownSampleFactor_xy_z(1),'bilinear'),hdr.transverseRotation,'bilinear') + info.RescaleIntercept;
			else
				vol1(:,:,floor((si-1)/hdr.DownSampleFactor_xy_z(2))+1) = vol1(:,:,floor((si-1)/hdr.DownSampleFactor_xy_z(2))+1) + ...
					info.RescaleSlope*imrotate(imresize(double(temp(:,:,ti)),1/hdr.DownSampleFactor_xy_z(1),'bilinear'),hdr.transverseRotation,'bilinear') + info.RescaleIntercept;
			end
		catch % deals with a case from TWH where the SPECT header looked like a gated/dynamic scan but had only one frame
			if isequal(size(temp),size(vol1))
				DCSlogmsg('  - * Warning! Image dimensions do not match header - completing as static scan. *');
				vol1 = info.RescaleSlope*imrotate(double(temp),hdr.transverseRotation,'bilinear') + info.RescaleIntercept;
			else
				error('Image dimensions are inconsistent with header')
			end
		end
		
	end % slice loop
	
	% Dead Time and Scatter
	if isnan(hdr.PETTotalCounts(fi)) && isfield(info,'Private_0009_10a7')
		hdr.PETTotalCounts(fi) = info.Private_0009_10a7;
	end
	if isnan(hdr.PrimaryPromptsCountsAccumulated(fi)) && isfield(info,'PrimaryPromptsCountsAccumulated')
		hdr.PrimaryPromptsCountsAccumulated(fi) = info.PrimaryPromptsCountsAccumulated;
	end
	if isnan(hdr.ScatterFractionFactor(fi)) && isfield(info,'ScatterFractionFactor')
		hdr.ScatterFractionFactor(fi) = info.ScatterFractionFactor;
	end
	if isnan(hdr.DeadTimeFactor(fi)) && isfield(info,'DeadTimeFactor')
		hdr.DeadTimeFactor(fi) = info.DeadTimeFactor;
	end
	
	% Frame start time duration
	if isnan(hdr.frame_len(fi))
		if strcmpi(hdr.image_type,'dynamic')
			if isfield(info,'ActualFrameDuration') && ~isempty(info.ActualFrameDuration)
				hdr.frame_len(fi) = info.ActualFrameDuration;
			elseif isfield(info,'RotationInformationSequence') && isfield(info.RotationInformationSequence,'Item_1') && isfield(info.RotationInformationSequence.Item_1,'ActualFrameDuration') && ~isempty(info.RotationInformationSequence.Item_1.ActualFrameDuration)
				hdr.frame_len(fi) = info.RotationInformationSequence.Item_1.ActualFrameDuration;
			end
		elseif strcmpi(hdr.image_type,'gated') && isfield(info,'FrameTime') && ~isempty(info.FrameTime)
			hdr.frame_len(fi) = info.FrameTime;
		end
	end
	
	% Quantification
	[vol(:,:,:,fi), hdr.quant_dynamic(fi)] = double2twobyte(vol1);
end % frame loop
DCSlogmsg('Done',1);

switch chklist.frame.mode
	case {'Explicit'}
		hdr.frame_start = reshape(chklist.frame.vals, hdr.nframes,1);
		if any(strcmpi(chklist.frame.field,{'AcquisitionTime','AcquisitionDate','ContentTime','ContentDate'}))
			hdr.frame_start = round((hdr.frame_start - hdr.frame_start(1))*24*60*60*1000); % convert to ms
			if hdr.nframes>1
				hdr.image_type = 'dynamic';
			end
		end
	case {'Index','Explicit Slot','Assumed'}
		if strcmpi(hdr.image_type,'gated')
			hdr.frame_start = (0:1/hdr.nframes:1-1/hdr.nframes)';
			if any(isnan(hdr.frame_len))
				hdr.frame_len = ones(hdr.nframes,1)/hdr.nframes;
			end
		else
			hdr.frame_start = (1:hdr.nframes)';
		end
end

% What if the frame length could not be populated (AMI did not
% follow the DICOM format - might be gaps between time-frames in CT/MR)
enableTFAlignmentTest = hdr.nframes>1;
if strcmpi(hdr.image_type,'dynamic') || strcmpi(hdr.image_type,'gated')
	if hdr.nframes == 1 % mislablelled study type
		hdr.image_type = 'static';
	else
		for i=1:hdr.nframes
			if isnan(hdr.frame_len(i))
				if strcmpi(hdr.modality,'CT') && isfield(info,'RevolutionTime')
					DCSlogmsg(['Frame length for frame ' num2str(i) ' estimated from CT revolution time.']);
					hdr.frame_len(i) = info.RevolutionTime*1000;
					enableTFAlignmentTest = false;
				else
					DCSlogmsg(['Frame length for frame ' num2str(i) ' could not be retreived from DICOM. Making automated estimate.']);
					if any(strcmpi(infodcm.Manufacturer, {'Advanced Molecular Imaging Inc.','AMI'}))
						if i==1
							hdr.frame_len(i) = hdr.frame_start(i)*2;
						else
							hdr.frame_len(i) = (hdr.frame_start(i) - (hdr.frame_start(i-1)+hdr.frame_len(i-1)/2))*2;
						end
					else
						if i~=hdr.nframes
							hdr.frame_len(i) = hdr.frame_start(i+1)-hdr.frame_start(i);
						else
							hdr.frame_len(i) = hdr.frame_len(i-1);
						end
					end
				end
			end
		end
	end
end

% Another interesting exception - Philips TF ActualFrameDuration field is not
% populated correctly
if enableTFAlignmentTest
	endTime = hdr.frame_start+hdr.frame_len;
	if ~all(hdr.frame_start(2:end)==endTime(1:end-1))
		DCSlogmsg('!!! Warning: Frame start and end times don''t line up. Infering frame lengths. !!!');
		hdr.frame_len(1:end-1) = hdr.frame_start(2:end) - hdr.frame_start(1:end-1);
		hdr.frame_len(end) = hdr.frame_len(end-1); % assume that last two frames are the same length.
	end
end


% Adjust image dimensions to down sampling
if any(hdr.DownSampleFactor_xy_z>1)
	% scale each pixel to the number of elements that it averages
	vol = vol/hdr.DownSampleFactor_xy_z(1)^2/hdr.DownSampleFactor_xy_z(2);
	if rem(hdr.xdim,hdr.DownSampleFactor_xy_z(1))>0
		vol(end,:,:,:) = vol(end,:,:,:)*hdr.DownSampleFactor_xy_z(1)/rem(hdr.xdim,hdr.DownSampleFactor_xy_z(1));
	end
	if rem(hdr.ydim,hdr.DownSampleFactor_xy_z(1))>0
		vol(:,end,:,:) = vol(:,end,:,:)*hdr.DownSampleFactor_xy_z(1)/rem(hdr.ydim,hdr.DownSampleFactor_xy_z(1));
	end
	if rem(hdr.nplanes,hdr.DownSampleFactor_xy_z(1))>0
		vol(:,:,end,:) = vol(:,:,end,:)*hdr.DownSampleFactor_xy_z(2)/rem(hdr.nplanes,hdr.DownSampleFactor_xy_z(2));
	end
	% adjust image dimensions and pixel sizes
	hdr.xdim = size(vol,1);
	hdr.ydim = size(vol,2);
	hdr.nplanes = size(vol,3);
	hdr.pix_mm_xy = hdr.pix_mm_xy*hdr.DownSampleFactor_xy_z(1);
	hdr.pix_mm_z = hdr.pix_mm_z*hdr.DownSampleFactor_xy_z(2);
end

%% Change image units
if isfield(infodcm,'Units')
	if strcmpi(infodcm.Units,'CNTS') && ~any(isnan(hdr.frame_len))
		if isfield(infodcm,'Private_7053_1009')
			hdr.quant_dynamic = hdr.quant_dynamic*info.Private_7053_1009; % convert to counts/sec (Phillips systems)
			hdr.image_units = 'Bq/cc';
		else
			hdr.quant_dynamic = hdr.quant_dynamic./(hdr.frame_len/1000); % convert to counts/sec (Phillips systems)
			hdr.image_units = 'CPS';
		end
	elseif strcmpi(infodcm.Units,'PROPCNTS') && ~any(isnan(hdr.frame_len))
		hdr.quant_dynamic = hdr.quant_dynamic./(hdr.frame_len/1000); % convert to counts/sec (Phillips systems)
		hdr.image_units = 'PROPCPS';
	end
elseif isfield(infodcm,'RealWorldValueMappingSequence') % found in Mediso NanoSPECT
		knownUnits = {'BECQUEREL','COUNTS'};
		fnames = fieldnames(infodcm.RealWorldValueMappingSequence);
		done = false;
		ki = 1;
		while ki<=length(knownUnits) && ~done
			for ti = 1:length(fnames)
				if isfield(infodcm.RealWorldValueMappingSequence.(fnames{ti}),'LUTLabel') &&...
						strcmpi(infodcm.RealWorldValueMappingSequence.(fnames{ti}).LUTLabel, knownUnits{ki})
					switch knownUnits{ki}
						case 'COUNTS'
							hdr.quant_dynamic = hdr.quant_dynamic * infodcm.RealWorldValueMappingSequence.(fnames{ti}).RealWorldValueSlope ...
								/ hdr.pix_mm_xy^2 / hdr.pix_mm_z * 1000; % per cc
							if any(isnan(hdr.frame_len))
								hdr.image_units = 'CPS';
							else
								hdr.quant_dynamic = hdr.quant_dynamic ./ (hdr.frame_len/1000);
								hdr.image_units = 'Bq/cc';
							end
							done = true;
						case 'BECQUEREL'
							hdr.quant_dynamic = hdr.quant_dynamic * infodcm.RealWorldValueMappingSequence.(fnames{ti}).RealWorldValueSlope...
								/ hdr.pix_mm_xy^2 / hdr.pix_mm_z * 1000; % per cc
							hdr.image_units = 'Bq/cc';
							done = true;
					end
				end
			end
		end
else
	if isfield(infodcm,'Manufacturer') && strcmpi(infodcm.Manufacturer,'GE MEDICAL SYSTEMS')
        if isfield(infodcm,'RescaleType') && strcmpi(infodcm.RescaleType,'HU')
			hdr.image_units = 'HU';
        else
            hdr.quant_dynamic = hdr.quant_dynamic ./ hdr.frame_len/1000 / (hdr.pix_mm_xy^2*hdr.pix_mm_z /1000);
            hdr.image_units = 'Bq/cc';
        end
	end
end

%% Accomodate Perfusion CT data from Siemens software
if isfield(info,'ImageComments')
	if strncmpi(info.ImageComments,'Blood Flow',10)
		DCSlogmsg('Detected a Siemens Blood Flow image - smoothing to remove NaN values')
		vol(vol<0) = nan;
		if ~isempty(strfind(info.ImageComments,'ml/100ml/min'))
			hdr.quant_dynamic = hdr.quant_dynamic/100; % Convert to ml/min/g units
		end
		hdr.image_units = 'mL/min/g';
		hdr.longitudinalFlip = false;
	end
end

%% Flip image as required to achieve standard anatomical orientation with head at slice 1
if ~hdr.longitudinalFlip 
	for i=1:hdr.nframes % one frame at a time to accomodate low memory.
		vol(:,:,:,i) = flip(vol(:,:,:,i),3);
	end
end