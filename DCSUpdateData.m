% DCSUpdateData - Update the series dicom data structure for the 
% DICOMConvertServer based the next image file.
%
% Usage:
% [DCSCurrentData valid] = DCSUpdateData(DCSCurrentData, source) - source 
% is one of the DICOM series files.
% [DCSCurrentData valid] = DCSUpdateData(DCSCurrentData, infodcm) - infodcm
% is one of the DICOM header structure.
%
% See also: dicomConvertServer, DCSInitData, DCSFinalizeData,
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


function [DCSCurrentData, valid] = DCSUpdateData(DCSCurrentData, infodcm)

if isstruct(infodcm)
	source =infodcm.Filename;
else
	source = infodcm;
	infodcm = dicominfo(source);
end

% Update the DCSCurrentData.chklist
switch DCSCurrentData.chklist.slice.mode
	case 'Explicit Location'
		slice = infodcm.(DCSCurrentData.chklist.slice.field);
	case 'Explicit Slices'
		try
			slice = DCSCurrentData.chklist.slice.fieldFactor*double(infodcm.(DCSCurrentData.chklist.slice.field));
		catch
			slice = DCSCurrentData.chklist.slice.fieldFactor*double(DCSCurrentData.infodcm.(DCSCurrentData.chklist.slice.field));
		end
	case 'Indexed'
		indx = infodcm.(DCSCurrentData.chklist.slice.field);
		if indx==0
			indx = 1;
		end
		if ~isnan(DCSCurrentData.hdr.nplanes)
			slice=mod(indx-1,DCSCurrentData.hdr.nplanes)+1;
			frame = floor(double(indx-1)/double(DCSCurrentData.hdr.nplanes))+1;
		else
			slice = floor((indx-1)/DCSCurrentData.hdr.nframes)+1;
			frame = mod(indx,DCSCurrentData.hdr.nframes)+1;
		end
	case 'Assumed', slice = 1;
end

switch DCSCurrentData.chklist.frame.mode
	case 'Explicit'
		if strcmpi(DCSCurrentData.chklist.frame.field,'ImageID')
			frame = getFrameNumFromID(infodcm.ImageID);
		elseif any(strcmpi(DCSCurrentData.chklist.frame.field,{'AcquisitionDate','ContentDate'}))
			frame = datenum(['00000101' infodcm.(strrep(DCSCurrentData.chklist.frame.field,'Date','Time'))],'yyyymmddHHMMSS') + ...
				(datenum(infodcm.(DCSCurrentData.chklist.frame.field),'yyyymmdd')-DCSCurrentData.chklist.frame.dateref);
		elseif strcmpi(DCSCurrentData.chklist.frame.field,'AcquisitionTime') || strcmpi(DCSCurrentData.chklist.frame.field,'ContentTime')
			frame = datenum(infodcm.(DCSCurrentData.chklist.frame.field),'HHMMSS');
		else
			frame = infodcm.(DCSCurrentData.chklist.frame.field);
			if ischar(frame)
				frame = str2double(frame);
			end
		end
	case 'Explicit Slot', frame = infodcm.(DCSCurrentData.chklist.frame.field);
	case 'Indexed' % already resolved in the slice section
	case 'Assumed', frame = 1;
end

% frame and slice are now resolved - update the checklist table
if ~isequal(size(frame),size(slice))
	if numel(frame)==1
		frame = frame*ones(size(slice));
	elseif numel(slice)==1
		slice = slice*ones(size(frame));
	end
end

valid = true;
for i=1:length(slice)
	si = getArrIndx(DCSCurrentData.chklist.slice.vals,slice(i));
	if si<=length(DCSCurrentData.chklist.slice.vals)
		if DCSCurrentData.chklist.slice.vals(si)~=slice(i)% a new plane
			DCSCurrentData.chklist.slice.vals = [DCSCurrentData.chklist.slice.vals(1:si-1) slice(i) DCSCurrentData.chklist.slice.vals(si:end)];
			DCSCurrentData.chklist.table = [DCSCurrentData.chklist.table(1:si-1,:,:); zeros(1,size(DCSCurrentData.chklist.table,2),2); DCSCurrentData.chklist.table(si:end,:,:)];
		end
	else
		DCSCurrentData.chklist.slice.vals = [DCSCurrentData.chklist.slice.vals slice(i)];
		DCSCurrentData.chklist.table = [DCSCurrentData.chklist.table; zeros(1,size(DCSCurrentData.chklist.table,2),2)];
	end
	
	fi = getArrIndx(DCSCurrentData.chklist.frame.vals,frame(i));
	if fi<=length(DCSCurrentData.chklist.frame.vals)
		if DCSCurrentData.chklist.frame.vals(fi)~=frame(i) % a new frame
			DCSCurrentData.chklist.frame.vals = [DCSCurrentData.chklist.frame.vals(1:fi-1) frame(i) DCSCurrentData.chklist.frame.vals(fi:end)];
			DCSCurrentData.chklist.table = [DCSCurrentData.chklist.table(:,1:fi-1,:), zeros(size(DCSCurrentData.chklist.table,1),1,2), DCSCurrentData.chklist.table(:,fi:end,:)];
		end
	else
		DCSCurrentData.chklist.frame.vals = [DCSCurrentData.chklist.frame.vals frame(i)];
		DCSCurrentData.chklist.table = [DCSCurrentData.chklist.table, zeros(size(DCSCurrentData.chklist.table,1),1,2)];
	end
	
	% No overlap detected - file valid?
	if DCSCurrentData.chklist.table(si,fi,1)
		DCSlogmsg([' *** overwrite conflict detected at ' source ' plane: ' num2str(si) ' frame: ' num2str(fi)]);
		valid = false;
	end
	
end




% add file to list and move to temporary directory
if valid
	DCSCurrentData.chklist.fileCount = DCSCurrentData.chklist.fileCount+1;
	for i=1:length(slice)
		fi = getArrIndx(DCSCurrentData.chklist.frame.vals,frame(i));
		si = getArrIndx(DCSCurrentData.chklist.slice.vals,slice(i));
		DCSCurrentData.chklist.table(si,fi,1) = DCSCurrentData.chklist.fileCount;
		DCSCurrentData.chklist.table(si,fi,2) = i;
	end
	if DCSCurrentData.chklist.useTempFiles
		DCSCurrentData.chklist.files = [DCSCurrentData.chklist.files source];
	else
		DCSCurrentData.chklist.files = [DCSCurrentData.chklist.files infodcm];
	end
end



%% For debugging purposese
% Checklist Table live display
if 0
	chklistFig = 102;
	if ~ishandle(chklistFig)
		figure(chklistFig);
		ah = axes;
		ih=imagesc(logical(DCSCurrentData.chklist.table(:,:,1)));
		xlabel(ah,'Frames'); ylabel(ah,'Slices');
		setappdata(chklistFig,'ImageHandle',ih);
		set(ah,'clim',[0 1]);
	else
		ih=getappdata(chklistFig,'ImageHandle');
		set(ih,'cdata',logical(DCSCurrentData.chklist.table(:,:,1)));
		set(get(ih,'parent'),'xlim',[0.5 size(DCSCurrentData.chklist.table,2)+0.5],'ylim',[0.5 size(DCSCurrentData.chklist.table,1)+0.5])
	end
	
end







%% HELPER FUNCTIONS

function i = getArrIndx(arr,val)
i = 1;
while i<=length(arr)
	if val>arr(i)
		i = i+1;
	else
		break;
	end
end

function frame = getFrameNumFromID(IDstring)
i = strfind(lower(IDstring),'frame');
IDstring = IDstring(i:end);
[~, b] = strtok(IDstring,' ');
b = strtok(b,' ');
frame = str2num(b);
