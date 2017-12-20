<<<<<<< HEAD
% Initialize the DICOM series internal data structure for the DICOMConvertServer
%
% Usage:
% DCSCurrentData = DCSInitCurrentData(source) - source is one of the DICOM
% series files.
%
% See also: dicomConvertServer, DCSUpdateData, DCSFinalizeData,
%           load_DICOMDirectory_scan, hdrInitDcm

% By Ran Klein August-2007
% University of Ottawa Heart Institute
% 2012-11-16 - RK - Partitioned this function out of dicomCenvertServer to
%                   be reused in load_DICOMDirectory_scan.


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function DCSCurrentData = DCSInitCurrentData(source)

[DCSCurrentData.hdr, DCSCurrentData.infodcm, messages] = hdrInitDcm(source); % get header information
for i=1:length(messages)
	DCSlogmsg(messages{i});
end
DCSCurrentData.chklist = [];

%% Determine the method for assigining slice (3rd dimension) in dynamic volume
if isfield(DCSCurrentData.infodcm,'ImplementationVersionName') && strcmpi(DCSCurrentData.infodcm.ImplementationVersionName, 'HERMESGOLD440')
	% Special case for Hermes download file as some genius overwrote the SliceLocation field with garbage
	DCSCurrentData.chklist.slice.mode = 'Indexed'; DCSCurrentData.chklist.slice.field = 'ImageIndex';
else
	if isfield(DCSCurrentData.infodcm,'SliceLocation')
		DCSCurrentData.chklist.slice.mode = 'Explicit Location'; DCSCurrentData.chklist.slice.field = 'SliceLocation';
		% 		elseif strcmpi(DCSCurrentData.infodcm.Manufacturer, 'GE MEDICAL SYSTEMS') && isfield(DCSCurrentData.infodcm,'Private_0009_1067')
		% 			DCSCurrentData.chklist.slice.mode = 'Explicit Location'; DCSCurrentData.chklist.slice.field = 'Private_0009_1067';
	elseif isfield(DCSCurrentData.infodcm,'SliceVector')
		DCSCurrentData.chklist.slice.mode = 'Explicit Slices'; DCSCurrentData.chklist.slice.field = 'SliceVector';
		if isfield(DCSCurrentData.infodcm,'SpacingBetweenSlices')
			DCSCurrentData.chklist.slice.fieldFactor = DCSCurrentData.infodcm.SpacingBetweenSlices;
		else
			DCSCurrentData.chklist.slice.fieldFactor = 1;
		end
	elseif ~isnan(DCSCurrentData.hdr.nplanes) && isfield(DCSCurrentData.infodcm,'ImageIndex')  % For Siemens IRW exported DICOM
		DCSCurrentData.chklist.slice.mode = 'Indexed'; DCSCurrentData.chklist.slice.field = 'ImageIndex';
	elseif ~isnan(DCSCurrentData.hdr.nplanes) && isfield(DCSCurrentData.infodcm,'InstanceNumber')  % For Siemens IRW exported DICOM
		DCSCurrentData.chklist.slice.mode = 'Indexed'; DCSCurrentData.chklist.slice.field = 'InstanceNumber';
	else
		DCSCurrentData.chklist.slice.mode = 'Assumed';
	end
end

%% Determine the method for assigining time frame (4th dimension) in dynamic volume
if contains(upper(DCSCurrentData.infodcm.Manufacturer),'PHILIPS') && strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'AcquisitionTime') && isfield(DCSCurrentData.infodcm,'AcquisitionDate')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionDate';
elseif isfield(DCSCurrentData.infodcm,'AcquisitionTime') && isfield(DCSCurrentData.infodcm,'AcquisitionDate')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionDate';
elseif isfield(DCSCurrentData.infodcm,'ContentTime') && isfield(DCSCurrentData.infodcm,'ContentDate')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'ContentDate';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'AcquisitionTime')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'ContentTime';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'ContentTime')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionTime';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'FrameReferenceTime') % Mid-frame based on decay
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'FrameReferenceTime';
elseif strcmpi(DCSCurrentData.hdr.image_type,'gated') && isfield(DCSCurrentData.infodcm,'TriggerTime')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'TriggerTime';
elseif isfield(DCSCurrentData.infodcm,'TimeSlotVector')
	DCSCurrentData.chklist.frame.mode = 'Explicit Slot'; DCSCurrentData.chklist.frame.field = 'TimeSlotVector';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'TemporalPositionIndex') && isfield(DCSCurrentData.infodcm,'TemporalPositionIndex')
	DCSCurrentData.chklist.frame.mode = 'Explicit Slot'; DCSCurrentData.chklist.frame.field = 'TemporalPositionIndex';
elseif strcmp(DCSCurrentData.chklist.slice.mode,'Indexed') && ~isnan(DCSCurrentData.hdr.nplanes)
	DCSCurrentData.chklist.frame.mode = 'Indexed';
elseif strcmpi(DCSCurrentData.hdr.image_type,'static') 
	if contains(upper(DCSCurrentData.infodcm.Manufacturer),'GE MEDICAL SYSTEMS') && isfield(DCSCurrentData.infodcm,'PctRpeakDelay')
		% Cardiac gated CT sequence of separately reconstructed images
		% (GE VCT) % added by RK - 2014-05-01
		DCSCurrentData.hdr.image_type = 'gated';
		DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'PctRpeakDelay';
	elseif contains(upper(DCSCurrentData.infodcm.Manufacturer),'GE MEDICAL SYSTEMS') && isfield(DCSCurrentData.infodcm,'AcquisitionTime')
		% CZT dynamic sequence of separately reconstructed images
		% (Exceleris 3.0562)
		DCSCurrentData.hdr.image_type = 'dynamic';
		DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionTime';
	elseif isfield(DCSCurrentData.infodcm,'ImageID')
		if ~isnan(getFrameNumFromID(DCSCurrentData.infodcm.ImageID)) % bioscan gated file???
			DCSCurrentData.hdr.image_type = 'gated';
			DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'ImageID';
		else
			DCSCurrentData.chklist.frame.mode = 'Assumed';
		end
	else
		DCSCurrentData.chklist.frame.mode = 'Assumed';
	end
else
	DCSCurrentData.chklist.frame.mode = 'Assumed';
end

%% Inititalize the checklist data
switch DCSCurrentData.chklist.slice.mode
	case 'Explicit Location'
		DCSCurrentData.chklist.slice.vals = sort(DCSCurrentData.infodcm.(DCSCurrentData.chklist.slice.field)); 
	case 'Explicit Slices'
		DCSCurrentData.chklist.slice.vals = sort(DCSCurrentData.chklist.slice.fieldFactor*double(unique(DCSCurrentData.infodcm.(DCSCurrentData.chklist.slice.field))));
	case 'Indexed'
		if isnan(DCSCurrentData.hdr.nplanes)
			DCSCurrentData.chklist.slice.vals = [];
		else
			DCSCurrentData.chklist.slice.vals = 1:DCSCurrentData.hdr.nplanes;
		end
	case 'Assumed'
		DCSCurrentData.chklist.slice.vals = 1;
end
switch DCSCurrentData.chklist.frame.mode
	case 'Explicit'
		if strcmpi(DCSCurrentData.chklist.frame.field,'ImageID')
			DCSCurrentData.chklist.frame.vals = getFrameNumFromID(DCSCurrentData.infodcm.ImageID);
		elseif any(strcmpi(DCSCurrentData.chklist.frame.field,{'AcquisitionDate','ContentDate'}))
			DCSCurrentData.chklist.frame.vals = datenum(['00000101' DCSCurrentData.infodcm.(strrep(DCSCurrentData.chklist.frame.field,'Date','Time'))],'yyyymmddHHMMSS');
			DCSCurrentData.chklist.frame.dateref = datenum(DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field),'yyyymmdd');
		elseif strcmpi(DCSCurrentData.chklist.frame.field,'AcquisitionTime') || strcmpi(DCSCurrentData.chklist.frame.field,'ContentTime')
			DCSCurrentData.chklist.frame.vals = datenum(DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field),'HHMMSS');
		else
			DCSCurrentData.chklist.frame.vals = DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field);
			if ischar(DCSCurrentData.chklist.frame.vals)
				DCSCurrentData.chklist.frame.vals = str2double(DCSCurrentData.chklist.frame.vals);
			end
		end
	case 'Explicit Slot'
		DCSCurrentData.chklist.frame.vals = unique(DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field));
	case 'Indexed'
		DCSCurrentData.chklist.frame.vals = 1:DCSCurrentData.hdr.nframes;
	case 'Assumed'
		DCSCurrentData.chklist.frame.vals = 1;
end
DCSCurrentData.chklist.table = zeros(length(DCSCurrentData.chklist.slice.vals),length(DCSCurrentData.chklist.frame.vals),2);
DCSCurrentData.chklist.fileCount = 0;
DCSCurrentData.chklist.files = {};




%% HELPER FUNCTIONS
function frame = getFrameNumFromID(IDstring)
i = strfind(lower(IDstring),'frame');
IDstring = IDstring(i:end);
[~, b] = strtok(IDstring,' ');
b = strtok(b,' ');
=======
% Initialize the DICOM series internal data structure for the DICOMConvertServer
%
% Usage:
% DCSCurrentData = DCSInitCurrentData(source) - source is one of the DICOM
% series files.
%
% See also: dicomConvertServer, DCSUpdateData, DCSFinalizeData,
%           load_DICOMDirectory_scan, hdrInitDcm

% By Ran Klein August-2007
% University of Ottawa Heart Institute
% 2012-11-16 - RK - Partitioned this function out of dicomCenvertServer to
%                   be reused in load_DICOMDirectory_scan.


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function DCSCurrentData = DCSInitCurrentData(source)

[DCSCurrentData.hdr, DCSCurrentData.infodcm, messages] = hdrInitDcm(source); % get header information
for i=1:length(messages)
	DCSlogmsg(messages{i});
end
DCSCurrentData.chklist = [];

%% Determine the method for assigining slice (3rd dimension) in dynamic volume
if isfield(DCSCurrentData.infodcm,'ImplementationVersionName') && strcmpi(DCSCurrentData.infodcm.ImplementationVersionName, 'HERMESGOLD440')
	% Special case for Hermes download file as some genius overwrote the SliceLocation field with garbage
	DCSCurrentData.chklist.slice.mode = 'Indexed'; DCSCurrentData.chklist.slice.field = 'ImageIndex';
else
	if isfield(DCSCurrentData.infodcm,'SliceLocation')
		DCSCurrentData.chklist.slice.mode = 'Explicit Location'; DCSCurrentData.chklist.slice.field = 'SliceLocation';
		% 		elseif strcmpi(DCSCurrentData.infodcm.Manufacturer, 'GE MEDICAL SYSTEMS') && isfield(DCSCurrentData.infodcm,'Private_0009_1067')
		% 			DCSCurrentData.chklist.slice.mode = 'Explicit Location'; DCSCurrentData.chklist.slice.field = 'Private_0009_1067';
	elseif isfield(DCSCurrentData.infodcm,'SliceVector')
		DCSCurrentData.chklist.slice.mode = 'Explicit Slices'; DCSCurrentData.chklist.slice.field = 'SliceVector';
		if isfield(DCSCurrentData.infodcm,'SpacingBetweenSlices')
			DCSCurrentData.chklist.slice.fieldFactor = DCSCurrentData.infodcm.SpacingBetweenSlices;
		else
			DCSCurrentData.chklist.slice.fieldFactor = 1;
		end
	elseif ~isnan(DCSCurrentData.hdr.nplanes) && isfield(DCSCurrentData.infodcm,'ImageIndex')  % For Siemens IRW exported DICOM
		DCSCurrentData.chklist.slice.mode = 'Indexed'; DCSCurrentData.chklist.slice.field = 'ImageIndex';
	elseif ~isnan(DCSCurrentData.hdr.nplanes) && isfield(DCSCurrentData.infodcm,'InstanceNumber')  % For Siemens IRW exported DICOM
		DCSCurrentData.chklist.slice.mode = 'Indexed'; DCSCurrentData.chklist.slice.field = 'InstanceNumber';
	else
		DCSCurrentData.chklist.slice.mode = 'Assumed';
	end
end

%% Determine the method for assigining time frame (4th dimension) in dynamic volume
if contains(upper(DCSCurrentData.infodcm.Manufacturer),'PHILIPS') && strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'AcquisitionTime') && isfield(DCSCurrentData.infodcm,'AcquisitionDate')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionDate';
elseif isfield(DCSCurrentData.infodcm,'AcquisitionTime') && isfield(DCSCurrentData.infodcm,'AcquisitionDate')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionDate';
elseif isfield(DCSCurrentData.infodcm,'ContentTime') && isfield(DCSCurrentData.infodcm,'ContentDate')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'ContentDate';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'AcquisitionTime')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'ContentTime';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'ContentTime')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionTime';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'FrameReferenceTime') % Mid-frame based on decay
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'FrameReferenceTime';
elseif strcmpi(DCSCurrentData.hdr.image_type,'gated') && isfield(DCSCurrentData.infodcm,'TriggerTime')
	DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'TriggerTime';
elseif isfield(DCSCurrentData.infodcm,'TimeSlotVector')
	DCSCurrentData.chklist.frame.mode = 'Explicit Slot'; DCSCurrentData.chklist.frame.field = 'TimeSlotVector';
elseif strcmpi(DCSCurrentData.hdr.image_type,'dynamic') && isfield(DCSCurrentData.infodcm,'TemporalPositionIndex') && isfield(DCSCurrentData.infodcm,'TemporalPositionIndex')
	DCSCurrentData.chklist.frame.mode = 'Explicit Slot'; DCSCurrentData.chklist.frame.field = 'TemporalPositionIndex';
elseif strcmp(DCSCurrentData.chklist.slice.mode,'Indexed') && ~isnan(DCSCurrentData.hdr.nplanes)
	DCSCurrentData.chklist.frame.mode = 'Indexed';
elseif strcmpi(DCSCurrentData.hdr.image_type,'static') 
	if contains(upper(DCSCurrentData.infodcm.Manufacturer),'GE MEDICAL SYSTEMS') && isfield(DCSCurrentData.infodcm,'PctRpeakDelay')
		% Cardiac gated CT sequence of separately reconstructed images
		% (GE VCT) % added by RK - 2014-05-01
		DCSCurrentData.hdr.image_type = 'gated';
		DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'PctRpeakDelay';
	elseif contains(upper(DCSCurrentData.infodcm.Manufacturer),'GE MEDICAL SYSTEMS') && isfield(DCSCurrentData.infodcm,'AcquisitionTime')
		% CZT dynamic sequence of separately reconstructed images
		% (Exceleris 3.0562)
		DCSCurrentData.hdr.image_type = 'dynamic';
		DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'AcquisitionTime';
	elseif isfield(DCSCurrentData.infodcm,'ImageID')
		if ~isnan(getFrameNumFromID(DCSCurrentData.infodcm.ImageID)) % bioscan gated file???
			DCSCurrentData.hdr.image_type = 'gated';
			DCSCurrentData.chklist.frame.mode = 'Explicit'; DCSCurrentData.chklist.frame.field = 'ImageID';
		else
			DCSCurrentData.chklist.frame.mode = 'Assumed';
		end
	else
		DCSCurrentData.chklist.frame.mode = 'Assumed';
	end
else
	DCSCurrentData.chklist.frame.mode = 'Assumed';
end

%% Inititalize the checklist data
switch DCSCurrentData.chklist.slice.mode
	case 'Explicit Location'
		DCSCurrentData.chklist.slice.vals = sort(DCSCurrentData.infodcm.(DCSCurrentData.chklist.slice.field)); 
	case 'Explicit Slices'
		DCSCurrentData.chklist.slice.vals = sort(DCSCurrentData.chklist.slice.fieldFactor*double(unique(DCSCurrentData.infodcm.(DCSCurrentData.chklist.slice.field))));
	case 'Indexed'
		if isnan(DCSCurrentData.hdr.nplanes)
			DCSCurrentData.chklist.slice.vals = [];
		else
			DCSCurrentData.chklist.slice.vals = 1:DCSCurrentData.hdr.nplanes;
		end
	case 'Assumed'
		DCSCurrentData.chklist.slice.vals = 1;
end
switch DCSCurrentData.chklist.frame.mode
	case 'Explicit'
		if strcmpi(DCSCurrentData.chklist.frame.field,'ImageID')
			DCSCurrentData.chklist.frame.vals = getFrameNumFromID(DCSCurrentData.infodcm.ImageID);
		elseif any(strcmpi(DCSCurrentData.chklist.frame.field,{'AcquisitionDate','ContentDate'}))
			DCSCurrentData.chklist.frame.vals = datenum(['00000101' DCSCurrentData.infodcm.(strrep(DCSCurrentData.chklist.frame.field,'Date','Time'))],'yyyymmddHHMMSS');
			DCSCurrentData.chklist.frame.dateref = datenum(DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field),'yyyymmdd');
		elseif strcmpi(DCSCurrentData.chklist.frame.field,'AcquisitionTime') || strcmpi(DCSCurrentData.chklist.frame.field,'ContentTime')
			DCSCurrentData.chklist.frame.vals = datenum(DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field),'HHMMSS');
		else
			DCSCurrentData.chklist.frame.vals = DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field);
			if ischar(DCSCurrentData.chklist.frame.vals)
				DCSCurrentData.chklist.frame.vals = str2double(DCSCurrentData.chklist.frame.vals);
			end
		end
	case 'Explicit Slot'
		DCSCurrentData.chklist.frame.vals = unique(DCSCurrentData.infodcm.(DCSCurrentData.chklist.frame.field));
	case 'Indexed'
		DCSCurrentData.chklist.frame.vals = 1:DCSCurrentData.hdr.nframes;
	case 'Assumed'
		DCSCurrentData.chklist.frame.vals = 1;
end
DCSCurrentData.chklist.table = zeros(length(DCSCurrentData.chklist.slice.vals),length(DCSCurrentData.chklist.frame.vals),2);
DCSCurrentData.chklist.fileCount = 0;
DCSCurrentData.chklist.files = {};




%% HELPER FUNCTIONS
function frame = getFrameNumFromID(IDstring)
i = strfind(lower(IDstring),'frame');
IDstring = IDstring(i:end);
[~, b] = strtok(IDstring,' ');
b = strtok(b,' ');
>>>>>>> f95da97b0f7b59512a505a4b1d7d27ebe583c1e6
frame = str2double(b);