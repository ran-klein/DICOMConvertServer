<<<<<<< HEAD
% dicomConvertServer - server that monitors a directory and converts dicom
% series into individual .mat files.
%
% Background: DICOM images can be awkward to handle due to differing 
% implementations and slow load times, especcially with large volume series
% consisting of hundred of files. To address these limitations for analysis
% in Matlab, the DICOM Convert Server automatically converts DICOM series
% into Matlab .mat.
% The server monitors an incoming directory for new DICOm files and then
% sorts them to derive a DICOM volume (can be 3D or 4D data). The converted
% files are stored in an output directory.
% The server runs as a background application and needs no user
% intervention accept for initial configuration using a configurations txt
% file. The server progress can be monitored with a live log display and
% also a log file.
%
% The converted .mat files contain the following variables:
% •	vol – a binary 4-dimensional image with dimension order x, y, z, and 
%         time. Each element encodes a single pixel intensity in 2 byte 
%         integer, that can be converted to native image unit by 
%         multiplying by a corresponding scaling factor in 
%         hdr.quant_dynamic. 
% •	hdr – in the image header and is a structure with the following fileds:
%     o	nativefile – the original image filename (before conversion to .mat)
%     o	filename – name of current filename (after conversion to .mat)
%     o	nframes – number of time frames (4th image dimension)
%     o	nplanes – number of image planes (3rd image dimension)
%     o	ydim – number of pixels in y-dimension (2nd image dimension)
%     o	xdim – number of pixels in x-dimension (1st image dimension)
%     o	pix_mm_xy – pixel size in mm in x- and y-dimensions
%     o	pix_mm_z' – pixel size in mm in z-dimension
%     o	resolution – estimated image resolution (based on reconstruction 
%           parameters) in mm
%     o	image_offset_mm – image position offset in x, y and z dimensions in
%           mm
%     o	reconstruction – text description of image reconstruction type
%     o	frame_start – time of start of each time frame (relative to 
%           beginning of scan) in ms
%     o	frame_len – length of each time frame in ms
%     o	quant_dynamic – scaling factor for each time frame image to convert 
%           from 2-byte integer format to floating point.
%     o	image_units – pixel intensity units (e.g. Bq/cc, HU)
%     o	PETTotalCounts – total number of PET counts in each time frame
%     o	PrimaryPromptsCountsAccumulated – number of accumulated prompt 
%           counts for each time frame
%     o	ScatterFractionFactor – scatter fraction factor for each time frame
%     o	DeadTimeFactor – system dead time factor for each time frame
%     o	patient_name – patient name
%     o	patientID – patient identification number
%     o	model_num – imaging system model
%     o	tracer – imaging tracer
%     o	study – study description
%     o	date – series date and time
%     o	patientDOB – patient date of birth
%     o	patientSex – patient sex
%     o	patientHeight – patient height in cm
%     o	patientWeight – patient weight in kg
%     o	injectedActivity – amount of tracer activity injected in MBq
%     o	studyID – study identifier number
%     o	examType – type of exam (e.g. rest, stress)
%     o	modality – imaging modality type (e.g. PT, NM, CT, MR)
%     o	image_type – image type (e.g. static, gated, dynamic)
%     o	transverseRotation – degrees rotation applied to the image data 
%           during conversion
%     o	longitudinalFlip – whether a longitudinal (z-dimension) flip was 
%           applied to the image during conversion
% •	dcm – DICOM header structure of the first file that was read (see 
%           hdr.Nativefile) during conversion.

%
% Usage:
% dicomConvertServer - start the server using default settings file.
% dicomConvertServer(flag) - specifies an operation for the server 
% 'start' - start the server (default if no missing).
% 'refresh' - scan the directory once.
% 'stop' - stop the server.
% 'help' - display DicomConvertServer help.
% dicomConvertServer(flag, optfile) - also specifies a settings file to
% load.
%
% See also: hdrInitDcm, DCSInitCurrentData, DCSUpdateData, DCSFinalizeData, 
%           DCSlogmsg

% By Ran Klein August-2007
% University of Ottawa Heart Institute


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function dicomConvertServer(flag, optfile)
global DCSVars  DCSCurrentData  

if nargin==0 || isempty(flag)
	flag = 'start';
end

if strcmpi(flag,'help')
	disp;
	return
end

% Start or refresh flags
if strcmpi(flag,'start') || strcmpi(flag,'refresh')
	if isfield(DCSVars,'logfileh') && ~isempty(DCSVars.logfileh)
		try fclose(DCSVars.logfileh); catch, end
	end
		
	%% Default options
	
	% the root directory based on the OS
	if strncmp(computer,'PC',2)
		rootdir = 'C:\';
	else  % isunix
		tdir = pwd;
		!cd $home
		roordir = pwd;
		cd(tdir);
	end
	DCSVars = struct('version','2.2',...
		'indir',[rootdir 'DICOMImport'],...
		'outdir',[rootdir 'DataStash'],...
		'logdir',[rootdir 'DataStash' filesep 'Log'],...
		'tempdir',[rootdir 'DataStash' filesep 'Temp'],...
		'cleanup',true,...
		'overwrite',false);

	%% Load options file if it exists
	if nargin<2
		if isdeployed
			path = pwd;
		else
			path = fileparts(mfilename('fullpath'));
		end
		optfile = [path filesep 'DICOMServer.dat'];
	end
	if exist(optfile,'file')
		of = fopen(optfile,'r');
		a = '';
		while ~isequal(a,-1)
			a = fgetl(of);
			if ~isequal(a,-1) && ~isempty(a) && a(1)~='%' % not a comment or empty line
				i = find(a=='='); 
				if length(i)==1 % valid lines only have one 
					temp = strtrim(a(i+1:end));
					if ~isempty(temp) % value is not empty
						if ~isempty(str2num(temp)) % value is a number
							eval(['DCSVars.' strtrim(a(1:i-1)) ' = ' temp ';']);
						else % value is a string
							eval(['DCSVars.' strtrim(a(1:i-1)) ' = ''' temp ''';']);
						end
					end % empty value
				end
			end % comment check
		end
		fclose(of);
	else
		optfile = '';
	end

	%% Startup log messages
	filename = [DCSVars.logdir filesep 'ConvertServerLog' datestr(now,'yyyy-mmm-dd') '.txt'];
	if ~exist(DCSVars.logdir,'dir')
		mkdir(DCSVars.logdir);
		mklogindir = true;
	else
		mklogindir = false;
	end
	DCSVars.logfileh = fopen(filename,'a');
	fprintf(DCSVars.logfileh,'\nStarting DICOM convert server session log: %s\n',filename);
	DCSlogmsg(['Starting DICOM Convert Server (V' DCSVars.version ')']);
	
	if mklogindir
		DCSlogmsg(['Logging directory created: ' DCSVars.logdir]);
	end
	
	if isempty(optfile)
		DCSlogmsg(['Warning: Opt file could not be found. ' optfile]);
		DCSlogmsg(' - Using default settings.');
	else
		DCSlogmsg(['Using options from opt file: ' optfile]);
	end
		
	DCSlogmsg('Settings:');
	DCSlogmsg('---------');
	DCSlogmsg(['Incoming directory: ' DCSVars.indir]);
	DCSlogmsg(['Temporary directory: ' DCSVars.tempdir]);
	DCSlogmsg(['Outgoing directory: ' DCSVars.outdir]);
	DCSlogmsg(['Logging directory: ' DCSVars.logdir]);
	if DCSVars.cleanup
		DCSlogmsg('Clean up incoming directory: yes');
	else
		DCSlogmsg('Clean up incoming directory: no');
	end
	if DCSVars.overwrite
		DCSlogmsg('Overwrite target files: yes');
	else
		DCSlogmsg('Overwrite target files: no');
	end
	DCSlogmsg('---------');
	
	% Check/Create missing directories
	if ~exist(DCSVars.indir,'dir')
		DCSlogmsg([' *** warning: Input directory ' DCSVars.indir ' was not found.']);
	end
	if ~exist(DCSVars.outdir,'dir')
		mkdir(DCSVars.outdir)
		DCSlogmsg(['Output directory created: ' DCSVars.outdir]);
	end
	if ~exist(DCSVars.tempdir,'dir')
		mkdir(DCSVars.tempdir)
		DCSlogmsg(['Temporary directory created: ' DCSVars.tempdir]);
	end

	%% Initialize data structures
	% dicom dictionary settings and respective codes
	DCSVars.dict = {'SeriesInstanceUID',...
		'ImageIndex',...
		'InstanceNumber',...
		'ImageID',...
		'RescaleSlope',...
		'SliceLocation',...
		'Private_0009_1067',... slice location on GE
		'RescaleIntercept',...
		'AcquisitionTime',...
		'AcquisitionDate',...
		'FrameReferenceTime',...
		'ActualFrameDuration',...
		'TriggerTime',... % gate start on GE
		'FrameTime',...%gate durations on GE
		'SliceVector',... % volumes in each file (SPECT)
		'TimeSlotVector',...
		'ScatterFractionFactor',...
		'DeadTimeFactor'};
	DCSVars.codes = uint32(zeros(1,length(DCSVars.dict)));
	for i=1:length(DCSVars.dict)
		[group, element] = dicomlookup(DCSVars.dict{i});
		if ~isempty(group)
			DCSVars.codes(i) = group + element*65536;
		else
			error(['Field ' DCSVars.dict{i} ' was not found in the global dictionary.']);
		end
	end
	
	DCSVars.currentfname = '';
	DCSVars.dirCheckCounter = 1;
	
	DCSCurrentData = struct('hdr',[],'infodcm',[],'chklist',[]);
	
	% Start flag - start the server as continuously monitoring
	if strcmpi(flag,'start')
		DCSVars.timerh = timer('TimerFcn',@directoryCheckWrapper,...
			'ExecutionMode','FixedSpacing','Period',5.0,'BusyMode','Drop');
		start(DCSVars.timerh);
	% refresh flag - single update cycle and done	
	else
		directoryCheckWrapper;
		try	fclose(DCSVars.logfileh);	catch, end
	end
	
	
	
% Stop flag	- the user closed the log window (loDCSlogmsg) indicating to
% terminate the server
elseif strcmpi(flag,'stop') 
	
	if isfield(DCSVars,'timerh') && ~isempty(DCSVars.timerh) && isvalid(DCSVars.timerh)
		stop(DCSVars.timerh);
		if strcmpi(get(DCSlogmsg,'UserData'),'Idle') % not in the middle of an operation.
			dicomConvertServer('Kill');
		end
	else % for some reason the timer could not be found
		DCSlogmsg('!!!! Timer could not be found')
		delete(DCSlogmsg);
	end
	
	
	
% Kill flag - clean up and terminate the server
elseif strcmpi(flag,'kill')
	
	savecurrentfname;
	try	delete(DCSVars.timerh);	catch, end
	delete(DCSlogmsg);
		
else
	warning(['Unknown directoryMonitor flag (' flag ') encountered']);
end



% Wrapper function to handle handshake operations for precessing status
function fcount = directoryCheckWrapper(varargin)
global DCSVars
set(DCSlogmsg,'UserData','Processing'); % Indicate that in middle of operation
fcount = directoryCheck(varargin);
if serverStopRequested(DCSVars)
	 dicomConvertServer('Kill');
else
	set(DCSlogmsg,'UserData','Idle'); % Indicate that not in middle of operation
end


% Routine for checking the incoming directory for new DICOMs to process
function fcount = directoryCheck(obj, event, string_arg)
global DCSVars 

fcount = 0; % have not processed any files in this directory
if nargin<1 || ~ischar(obj)
	DCSlogmsg('Checking incoming directory...');
	dir = DCSVars.indir;
	ignorelist = {'.','..','DICOMDIR'};
else
	dir = obj;
	ignorelist = {'.','..'};
end
files = listfiles('*.*',dir,'d');

% screen eachfile in the incoming directory
for i=1:length(files)
	if serverStopRequested(DCSVars)
		break;
	end
	if ~any(strcmpi(ignorelist,files{i}))
		if isdir([dir filesep files{i}])
			dirfcount = directoryCheck([dir filesep files{i}]);
			if serverStopRequested(DCSVars)
				break;
			end
			if dirfcount==0 && ... no file was processed
					DCSVars.cleanup && isempty(listfiles('*.*',[dir filesep files{i}],'d')) % has the entire directory been processed?
				DCSlogmsg(['Removing: ' dir filesep files{i}]);
				try
					rmdir([dir filesep files{i}]);
				catch
				end
			end
		else % not a directory - then a file
			try
				sortDicom([dir filesep files{i}]); % screen each file
				fcount = fcount+1; % another file was processed
			catch
				DCSlogmsg('  *** Error detected while sorting');
				if DCSVars.cleanup
					delete([dir filesep files{i}]);
				end
			end
		end
	end
end % files loop

% Purge and cleanup of temp files that are incomplete or suspected as
% complete in cases where the number of frames is unknown.
if ~serverStopRequested(DCSVars) && strcmp(dir,DCSVars.indir)
	savecurrentfname;
	DCSVars.currentfname = '';
	DCSVars.dirCheckCounter = DCSVars.dirCheckCounter+1;
	if DCSVars.dirCheckCounter>5
		DCSVars.dirCheckCounter=1;
	end
	files = listfiles('*.mat',DCSVars.tempdir);
	if ~serverStopRequested(DCSVars) && ~isempty(files)
		DCSlogmsg('Processing sorted files...');
		for i=1:length(files)
			if serverStopRequested(DCSVars)
				break;
			end
			if ~strcmpi(files{i},'.') && ~strcmpi(files{i},'..')
				t = load([DCSVars.tempdir filesep files{i}],'chklist','dirCheckCounter');
				if t.dirCheckCounter == DCSVars.dirCheckCounter % counter elapsed - deal with this file
					fname = [DCSVars.tempdir filesep files{i}(1:end-4)];
					if isfield(t.chklist,'table') && all(all(t.chklist.table(:,:,1))) % the file is complete
						DCSlogmsg(['Accepting ' files{i}]);
						try
							savetarget(fname);
							if serverStopRequested(DCSVars)
								break;
							end
						catch l
							DCSlogmsg('  *** Error while processing');
							DCSlogmsg(['      - ' l.message]);
							if strcmp(l.identifier,'MATLAB:nomem') % Added by RK - 2010-10-12 to deal with run out of memory error
								DCSlogmsg('      -    Preserving temp files for future attempt.');
								continue
							end
						end
					else
						DCSlogmsg(['Rejecting ' files{i}]);
						DCSlogmsg([' - Missing files: ' num2str(length(t.chklist.slice.vals)) ' slices X '...
							num2str(length(t.chklist.frame.vals)) ' frames = '...
							num2str(length(t.chklist.slice.vals) * length(t.chklist.frame.vals)) ' images, but only '...
							num2str(max(max(t.chklist.table(:,:,1))) * max(max(t.chklist.table(:,:,2)))) ' images found']);
					end
					delete([fname '.mat']);
					if exist(fname,'dir')
                        % Changed by RK 2011-03-01 from '*.*' to '*'
						delete([fname filesep '*']);
						rmdir(fname);
					end
				end
			end
		end % files loop
		DCSlogmsg('Processing cycle complete.');
	end % any temp files
end


%% Routine to sort the dicom files
function sortDicom(source)
global DCSVars  DCSCurrentData 
DCSlogmsg(['Sorting: ' source]);
try
	info = dicominfo(source);
	if strcmpi(info.Manufacturer, 'TOSHIBA') &&...
			isfield(info,'ScanOptions') && strcmpi(info.ScanOptions,'DVOLUME_CT')
		% addresses Toshiba dynamic CT using different series UID for each
		% time frame.
		info.SeriesInstanceUID = info.SeriesInstanceUID(1:find(info.SeriesInstanceUID=='.',1,'last')-1);
	end
	fname = ['temp_' strrep(info.SeriesInstanceUID,'.','_')];
catch
	if DCSVars.cleanup
		DCSlogmsg(' - Removing INVALID source file');
		% Note! - At the moment, files are deleted, if they are support
		% files they may need to be copied over/preserved/linked????
		try
			delete(source);
		catch
			DCSlogmsg('   - file removal failed');
		end
	end
	return
end

if ~strcmpi(DCSVars.currentfname,fname) % we've encountered a new image file
	savecurrentfname; % save the current file
	DCSVars.currentfname = fname;
	if exist([DCSVars.tempdir filesep DCSVars.currentfname '.mat'],'file') % has a temp file already been created?
		DCSlogmsg([' - Reverting to previous series: ' fname]);
		DCSCurrentData = load([DCSVars.tempdir filesep fname '.mat'],'-mat'); % load the new current file to append data
		DCSCurrentData = rmfield(DCSCurrentData,'dirCheckCounter'); % restore the structure format to compensate for savecurrentfname
	else % Create a new file
		DCSlogmsg(' - Detected new series');
		
		DCSCurrentData = DCSInitCurrentData(source);
	end
else
	DCSlogmsg(' - Adding to current series');
end

[DCSCurrentData, valid] = DCSUpdateData(DCSCurrentData, info);

% updae the filename in the list and move to temporary directory
if valid
	path = [DCSVars.tempdir filesep 'temp_' strrep(info.SeriesInstanceUID,'.','_')];
	if ~exist(path,'dir'),	mkdir(path); end
	tfname = sprintf('%d',100000+DCSCurrentData.chklist.fileCount);
	tfname = [path filesep 'IMG' tfname(2:end)];
	if DCSVars.cleanup
		movefile(source, tfname, 'f');
	else
		copyfile(source, tfname, 'f');
	end
	DCSCurrentData.chklist.files{end} = tfname;

else
	if DCSVars.cleanup
		delete(source);
	end
end


%% Save the current file being processed
function savecurrentfname
global DCSVars DCSCurrentData %#ok<NUSED>
if ~isempty(DCSVars.currentfname)
	try
		save([DCSVars.tempdir filesep DCSVars.currentfname '.mat'],'-mat','-struct','DCSCurrentData'); % Do not save as a structur so that can quickly read only ncessary fields
		dirCheckCounter = DCSVars.dirCheckCounter; %#ok<NASGU>
		save([DCSVars.tempdir filesep DCSVars.currentfname '.mat'],'-mat','-append','dirCheckCounter');
	catch
		DCSlogmsg(['Failed while saving: ' DCSVars.tempdir filesep DCSVars.currentfname '.mat']);
	end
end


%% Save the completely processed file to the target
function savetarget(source)
global DCSVars % Don't need global DCSCurrentData since comes from source file

load([source '.mat'],'hdr','infodcm','chklist');

[vol, hdr] = DCSFinalizeData(struct('infodcm',infodcm,'chklist',chklist,'hdr',hdr)); %#ok<ASGLU,NODEF>

if strcmpi(hdr.patient_name,hdr.patientID)
	fname = [titleCase(hdr.patient_name) '_' findTracer({hdr, 'UnknownTracer'}) '_' titleCase(hdr.study)];
else
	fname = [titleCase(hdr.patient_name) '_' hdr.patientID '_' findTracer({hdr, 'UnknownTracer'}) '_' titleCase(hdr.study)];
end

i = 1;
while i<=length(fname)
	if (fname(i)>='A' && fname(i)<='Z') || ...
			(fname(i)>='a' && fname(i)<='z') || ...
			(fname(i)>='0' && fname(i)<='9') || ...
			any(fname(i)=='_()@')
		i = i+1;
	else
		fname(i)='';
	end
end
fname = [DCSVars.outdir filesep fname];

if ~DCSVars.overwrite && exist([fname '.mat'],'file')
	i = 1;
	while exist([fname '_' num2str(i) '.mat'],'file')
		i = i+1;
	end
	fname = [fname '_' num2str(i) '.mat'];
else
	fname = [fname '.mat'];
end
		
hdr.DICOMConversionNotes = ['Converted ' datestr(now,'yyyy-mm-dd HH:MM:SS') ' by ' mfilename ' (V' DCSVars.version ')'];
try
	save(fname,'-mat','-v7.3','hdr','vol','infodcm');
	DCSlogmsg(['Saved to:' fname]);
catch
	DCSlogmsg(' - Error while saving');
	savecurrentfname
	return
end
DCSVars.currentfname = ''; % this is no longer the current file to process


%% Switch string to title case
function str = titleCase(str)
str = lower(str);
indx = [0 find(str==' ')];
for i=1:length(indx)
	if length(str)>indx(i)
		str(indx(i)+1) = upper(str(indx(i)+1));
	end
end


%% Checks if user requested to close the server
function status = serverStopRequested(DCSVars)
=======
% dicomConvertServer - server that monitors a directory and converts dicom
% series into individual .mat files.
%
% Background: DICOM images can be awkward to handle due to differing 
% implementations and slow load times, especcially with large volume series
% consisting of hundred of files. To address these limitations for analysis
% in Matlab, the DICOM Convert Server automatically converts DICOM series
% into Matlab .mat.
% The server monitors an incoming directory for new DICOm files and then
% sorts them to derive a DICOM volume (can be 3D or 4D data). The converted
% files are stored in an output directory.
% The server runs as a background application and needs no user
% intervention accept for initial configuration using a configurations txt
% file. The server progress can be monitored with a live log display and
% also a log file.
%
% The converted .mat files contain the following variables:
% •	vol – a binary 4-dimensional image with dimension order x, y, z, and 
%         time. Each element encodes a single pixel intensity in 2 byte 
%         integer, that can be converted to native image unit by 
%         multiplying by a corresponding scaling factor in 
%         hdr.quant_dynamic. 
% •	hdr – in the image header and is a structure with the following fileds:
%     o	nativefile – the original image filename (before conversion to .mat)
%     o	filename – name of current filename (after conversion to .mat)
%     o	nframes – number of time frames (4th image dimension)
%     o	nplanes – number of image planes (3rd image dimension)
%     o	ydim – number of pixels in y-dimension (2nd image dimension)
%     o	xdim – number of pixels in x-dimension (1st image dimension)
%     o	pix_mm_xy – pixel size in mm in x- and y-dimensions
%     o	pix_mm_z' – pixel size in mm in z-dimension
%     o	resolution – estimated image resolution (based on reconstruction 
%           parameters) in mm
%     o	image_offset_mm – image position offset in x, y and z dimensions in
%           mm
%     o	reconstruction – text description of image reconstruction type
%     o	frame_start – time of start of each time frame (relative to 
%           beginning of scan) in ms
%     o	frame_len – length of each time frame in ms
%     o	quant_dynamic – scaling factor for each time frame image to convert 
%           from 2-byte integer format to floating point.
%     o	image_units – pixel intensity units (e.g. Bq/cc, HU)
%     o	PETTotalCounts – total number of PET counts in each time frame
%     o	PrimaryPromptsCountsAccumulated – number of accumulated prompt 
%           counts for each time frame
%     o	ScatterFractionFactor – scatter fraction factor for each time frame
%     o	DeadTimeFactor – system dead time factor for each time frame
%     o	patient_name – patient name
%     o	patientID – patient identification number
%     o	model_num – imaging system model
%     o	tracer – imaging tracer
%     o	study – study description
%     o	date – series date and time
%     o	patientDOB – patient date of birth
%     o	patientSex – patient sex
%     o	patientHeight – patient height in cm
%     o	patientWeight – patient weight in kg
%     o	injectedActivity – amount of tracer activity injected in MBq
%     o	studyID – study identifier number
%     o	examType – type of exam (e.g. rest, stress)
%     o	modality – imaging modality type (e.g. PT, NM, CT, MR)
%     o	image_type – image type (e.g. static, gated, dynamic)
%     o	transverseRotation – degrees rotation applied to the image data 
%           during conversion
%     o	longitudinalFlip – whether a longitudinal (z-dimension) flip was 
%           applied to the image during conversion
% •	dcm – DICOM header structure of the first file that was read (see 
%           hdr.Nativefile) during conversion.

%
% Usage:
% dicomConvertServer - start the server using default settings file.
% dicomConvertServer(flag) - specifies an operation for the server 
% 'start' - start the server (default if no missing).
% 'refresh' - scan the directory once.
% 'stop' - stop the server.
% 'help' - display DicomConvertServer help.
% dicomConvertServer(flag, optfile) - also specifies a settings file to
% load.
%
% See also: hdrInitDcm, DCSInitCurrentData, DCSUpdateData, DCSFinalizeData, 
%           DCSlogmsg

% By Ran Klein August-2007
% University of Ottawa Heart Institute


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function dicomConvertServer(flag, optfile)
global DCSVars  DCSCurrentData  

if nargin==0 || isempty(flag)
	flag = 'start';
end

if strcmpi(flag,'help')
	disp;
	return
end

% Start or refresh flags
if strcmpi(flag,'start') || strcmpi(flag,'refresh')
	if isfield(DCSVars,'logfileh') && ~isempty(DCSVars.logfileh)
		try fclose(DCSVars.logfileh); catch, end
	end
		
	%% Default options
	
	% the root directory based on the OS
	if strncmp(computer,'PC',2)
		rootdir = 'C:\';
	else  % isunix
		tdir = pwd;
		!cd $home
		roordir = pwd;
		cd(tdir);
	end
	DCSVars = struct('version','2.2',...
		'indir',[rootdir 'DICOMImport'],...
		'outdir',[rootdir 'DataStash'],...
		'logdir',[rootdir 'DataStash' filesep 'Log'],...
		'tempdir',[rootdir 'DataStash' filesep 'Temp'],...
		'cleanup',true,...
		'overwrite',false);

	%% Load options file if it exists
	if nargin<2
		if isdeployed
			path = pwd;
		else
			path = fileparts(mfilename('fullpath'));
		end
		optfile = [path filesep 'DICOMServer.dat'];
	end
	if exist(optfile,'file')
		of = fopen(optfile,'r');
		a = '';
		while ~isequal(a,-1)
			a = fgetl(of);
			if ~isequal(a,-1) && ~isempty(a) && a(1)~='%' % not a comment or empty line
				i = find(a=='='); 
				if length(i)==1 % valid lines only have one 
					temp = strtrim(a(i+1:end));
					if ~isempty(temp) % value is not empty
						if ~isempty(str2num(temp)) % value is a number
							eval(['DCSVars.' strtrim(a(1:i-1)) ' = ' temp ';']);
						else % value is a string
							eval(['DCSVars.' strtrim(a(1:i-1)) ' = ''' temp ''';']);
						end
					end % empty value
				end
			end % comment check
		end
		fclose(of);
	else
		optfile = '';
	end

	%% Startup log messages
	filename = [DCSVars.logdir filesep 'ConvertServerLog' datestr(now,'yyyy-mmm-dd') '.txt'];
	if ~exist(DCSVars.logdir,'dir')
		mkdir(DCSVars.logdir);
		mklogindir = true;
	else
		mklogindir = false;
	end
	DCSVars.logfileh = fopen(filename,'a');
	fprintf(DCSVars.logfileh,'\nStarting DICOM convert server session log: %s\n',filename);
	DCSlogmsg(['Starting DICOM Convert Server (V' DCSVars.version ')']);
	
	if mklogindir
		DCSlogmsg(['Logging directory created: ' DCSVars.logdir]);
	end
	
	if isempty(optfile)
		DCSlogmsg(['Warning: Opt file could not be found. ' optfile]);
		DCSlogmsg(' - Using default settings.');
	else
		DCSlogmsg(['Using options from opt file: ' optfile]);
	end
		
	DCSlogmsg('Settings:');
	DCSlogmsg('---------');
	DCSlogmsg(['Incoming directory: ' DCSVars.indir]);
	DCSlogmsg(['Temporary directory: ' DCSVars.tempdir]);
	DCSlogmsg(['Outgoing directory: ' DCSVars.outdir]);
	DCSlogmsg(['Logging directory: ' DCSVars.logdir]);
	if DCSVars.cleanup
		DCSlogmsg('Clean up incoming directory: yes');
	else
		DCSlogmsg('Clean up incoming directory: no');
	end
	if DCSVars.overwrite
		DCSlogmsg('Overwrite target files: yes');
	else
		DCSlogmsg('Overwrite target files: no');
	end
	DCSlogmsg('---------');
	
	% Check/Create missing directories
	if ~exist(DCSVars.indir,'dir')
		DCSlogmsg([' *** warning: Input directory ' DCSVars.indir ' was not found.']);
	end
	if ~exist(DCSVars.outdir,'dir')
		mkdir(DCSVars.outdir)
		DCSlogmsg(['Output directory created: ' DCSVars.outdir]);
	end
	if ~exist(DCSVars.tempdir,'dir')
		mkdir(DCSVars.tempdir)
		DCSlogmsg(['Temporary directory created: ' DCSVars.tempdir]);
	end

	%% Initialize data structures
	% dicom dictionary settings and respective codes
	DCSVars.dict = {'SeriesInstanceUID',...
		'ImageIndex',...
		'InstanceNumber',...
		'ImageID',...
		'RescaleSlope',...
		'SliceLocation',...
		'Private_0009_1067',... slice location on GE
		'RescaleIntercept',...
		'AcquisitionTime',...
		'AcquisitionDate',...
		'FrameReferenceTime',...
		'ActualFrameDuration',...
		'TriggerTime',... % gate start on GE
		'FrameTime',...%gate durations on GE
		'SliceVector',... % volumes in each file (SPECT)
		'TimeSlotVector',...
		'ScatterFractionFactor',...
		'DeadTimeFactor'};
	DCSVars.codes = uint32(zeros(1,length(DCSVars.dict)));
	for i=1:length(DCSVars.dict)
		[group, element] = dicomlookup(DCSVars.dict{i});
		if ~isempty(group)
			DCSVars.codes(i) = group + element*65536;
		else
			error(['Field ' DCSVars.dict{i} ' was not found in the global dictionary.']);
		end
	end
	
	DCSVars.currentfname = '';
	DCSVars.dirCheckCounter = 1;
	
	DCSCurrentData = struct('hdr',[],'infodcm',[],'chklist',[]);
	
	% Start flag - start the server as continuously monitoring
	if strcmpi(flag,'start')
		DCSVars.timerh = timer('TimerFcn',@directoryCheckWrapper,...
			'ExecutionMode','FixedSpacing','Period',5.0,'BusyMode','Drop');
		start(DCSVars.timerh);
	% refresh flag - single update cycle and done	
	else
		directoryCheckWrapper;
		try	fclose(DCSVars.logfileh);	catch, end
	end
	
	
	
% Stop flag	- the user closed the log window (loDCSlogmsg) indicating to
% terminate the server
elseif strcmpi(flag,'stop') 
	
	if isfield(DCSVars,'timerh') && ~isempty(DCSVars.timerh) && isvalid(DCSVars.timerh)
		stop(DCSVars.timerh);
		if strcmpi(get(DCSlogmsg,'UserData'),'Idle') % not in the middle of an operation.
			dicomConvertServer('Kill');
		end
	else % for some reason the timer could not be found
		DCSlogmsg('!!!! Timer could not be found')
		delete(DCSlogmsg);
	end
	
	
	
% Kill flag - clean up and terminate the server
elseif strcmpi(flag,'kill')
	
	savecurrentfname;
	try	delete(DCSVars.timerh);	catch, end
	delete(DCSlogmsg);
		
else
	warning(['Unknown directoryMonitor flag (' flag ') encountered']);
end



% Wrapper function to handle handshake operations for precessing status
function fcount = directoryCheckWrapper(varargin)
global DCSVars
set(DCSlogmsg,'UserData','Processing'); % Indicate that in middle of operation
fcount = directoryCheck(varargin);
if serverStopRequested(DCSVars)
	 dicomConvertServer('Kill');
else
	set(DCSlogmsg,'UserData','Idle'); % Indicate that not in middle of operation
end


% Routine for checking the incoming directory for new DICOMs to process
function fcount = directoryCheck(obj, event, string_arg)
global DCSVars 

fcount = 0; % have not processed any files in this directory
if nargin<1 || ~ischar(obj)
	DCSlogmsg('Checking incoming directory...');
	dir = DCSVars.indir;
	ignorelist = {'.','..','DICOMDIR'};
else
	dir = obj;
	ignorelist = {'.','..'};
end
files = listfiles('*.*',dir,'d');

% screen eachfile in the incoming directory
for i=1:length(files)
	if serverStopRequested(DCSVars)
		break;
	end
	if ~any(strcmpi(ignorelist,files{i}))
		if isdir([dir filesep files{i}])
			dirfcount = directoryCheck([dir filesep files{i}]);
			if serverStopRequested(DCSVars)
				break;
			end
			if dirfcount==0 && ... no file was processed
					DCSVars.cleanup && isempty(listfiles('*.*',[dir filesep files{i}],'d')) % has the entire directory been processed?
				DCSlogmsg(['Removing: ' dir filesep files{i}]);
				try
					rmdir([dir filesep files{i}]);
				catch
				end
			end
		else % not a directory - then a file
			try
				sortDicom([dir filesep files{i}]); % screen each file
				fcount = fcount+1; % another file was processed
			catch
				DCSlogmsg('  *** Error detected while sorting');
				if DCSVars.cleanup
					delete([dir filesep files{i}]);
				end
			end
		end
	end
end % files loop

% Purge and cleanup of temp files that are incomplete or suspected as
% complete in cases where the number of frames is unknown.
if ~serverStopRequested(DCSVars) && strcmp(dir,DCSVars.indir)
	savecurrentfname;
	DCSVars.currentfname = '';
	DCSVars.dirCheckCounter = DCSVars.dirCheckCounter+1;
	if DCSVars.dirCheckCounter>5
		DCSVars.dirCheckCounter=1;
	end
	files = listfiles('*.mat',DCSVars.tempdir);
	if ~serverStopRequested(DCSVars) && ~isempty(files)
		DCSlogmsg('Processing sorted files...');
		for i=1:length(files)
			if serverStopRequested(DCSVars)
				break;
			end
			if ~strcmpi(files{i},'.') && ~strcmpi(files{i},'..')
				t = load([DCSVars.tempdir filesep files{i}],'chklist','dirCheckCounter');
				if t.dirCheckCounter == DCSVars.dirCheckCounter % counter elapsed - deal with this file
					fname = [DCSVars.tempdir filesep files{i}(1:end-4)];
					if isfield(t.chklist,'table') && all(all(t.chklist.table(:,:,1))) % the file is complete
						DCSlogmsg(['Accepting ' files{i}]);
						try
							savetarget(fname);
							if serverStopRequested(DCSVars)
								break;
							end
						catch l
							DCSlogmsg('  *** Error while processing');
							DCSlogmsg(['      - ' l.message]);
							if strcmp(l.identifier,'MATLAB:nomem') % Added by RK - 2010-10-12 to deal with run out of memory error
								DCSlogmsg('      -    Preserving temp files for future attempt.');
								continue
							end
						end
					else
						DCSlogmsg(['Rejecting ' files{i}]);
						DCSlogmsg([' - Missing files: ' num2str(length(t.chklist.slice.vals)) ' slices X '...
							num2str(length(t.chklist.frame.vals)) ' frames = '...
							num2str(length(t.chklist.slice.vals) * length(t.chklist.frame.vals)) ' images, but only '...
							num2str(max(max(t.chklist.table(:,:,1))) * max(max(t.chklist.table(:,:,2)))) ' images found']);
					end
					delete([fname '.mat']);
					if exist(fname,'dir')
                        % Changed by RK 2011-03-01 from '*.*' to '*'
						delete([fname filesep '*']);
						rmdir(fname);
					end
				end
			end
		end % files loop
		DCSlogmsg('Processing cycle complete.');
	end % any temp files
end


%% Routine to sort the dicom files
function sortDicom(source)
global DCSVars  DCSCurrentData 
DCSlogmsg(['Sorting: ' source]);
try
	info = dicominfo(source);
	if strcmpi(info.Manufacturer, 'TOSHIBA') &&...
			isfield(info,'ScanOptions') && strcmpi(info.ScanOptions,'DVOLUME_CT')
		% addresses Toshiba dynamic CT using different series UID for each
		% time frame.
		info.SeriesInstanceUID = info.SeriesInstanceUID(1:find(info.SeriesInstanceUID=='.',1,'last')-1);
	end
	fname = ['temp_' strrep(info.SeriesInstanceUID,'.','_')];
catch
	if DCSVars.cleanup
		DCSlogmsg(' - Removing INVALID source file');
		% Note! - At the moment, files are deleted, if they are support
		% files they may need to be copied over/preserved/linked????
		try
			delete(source);
		catch
			DCSlogmsg('   - file removal failed');
		end
	end
	return
end

if ~strcmpi(DCSVars.currentfname,fname) % we've encountered a new image file
	savecurrentfname; % save the current file
	DCSVars.currentfname = fname;
	if exist([DCSVars.tempdir filesep DCSVars.currentfname '.mat'],'file') % has a temp file already been created?
		DCSlogmsg([' - Reverting to previous series: ' fname]);
		DCSCurrentData = load([DCSVars.tempdir filesep fname '.mat'],'-mat'); % load the new current file to append data
		DCSCurrentData = rmfield(DCSCurrentData,'dirCheckCounter'); % restore the structure format to compensate for savecurrentfname
	else % Create a new file
		DCSlogmsg(' - Detected new series');
		
		DCSCurrentData = DCSInitCurrentData(source);
	end
else
	DCSlogmsg(' - Adding to current series');
end

[DCSCurrentData, valid] = DCSUpdateData(DCSCurrentData, info);

% updae the filename in the list and move to temporary directory
if valid
	path = [DCSVars.tempdir filesep 'temp_' strrep(info.SeriesInstanceUID,'.','_')];
	if ~exist(path,'dir'),	mkdir(path); end
	tfname = sprintf('%d',100000+DCSCurrentData.chklist.fileCount);
	tfname = [path filesep 'IMG' tfname(2:end)];
	if DCSVars.cleanup
		movefile(source, tfname, 'f');
	else
		copyfile(source, tfname, 'f');
	end
	DCSCurrentData.chklist.files{end} = tfname;

else
	if DCSVars.cleanup
		delete(source);
	end
end


%% Save the current file being processed
function savecurrentfname
global DCSVars DCSCurrentData %#ok<NUSED>
if ~isempty(DCSVars.currentfname)
	try
		save([DCSVars.tempdir filesep DCSVars.currentfname '.mat'],'-mat','-struct','DCSCurrentData'); % Do not save as a structur so that can quickly read only ncessary fields
		dirCheckCounter = DCSVars.dirCheckCounter; %#ok<NASGU>
		save([DCSVars.tempdir filesep DCSVars.currentfname '.mat'],'-mat','-append','dirCheckCounter');
	catch
		DCSlogmsg(['Failed while saving: ' DCSVars.tempdir filesep DCSVars.currentfname '.mat']);
	end
end


%% Save the completely processed file to the target
function savetarget(source)
global DCSVars % Don't need global DCSCurrentData since comes from source file

load([source '.mat'],'hdr','infodcm','chklist');

[vol, hdr] = DCSFinalizeData(struct('infodcm',infodcm,'chklist',chklist,'hdr',hdr)); %#ok<ASGLU,NODEF>

if strcmpi(hdr.patient_name,hdr.patientID)
	fname = [titleCase(hdr.patient_name) '_' findTracer({hdr, 'UnknownTracer'}) '_' titleCase(hdr.study)];
else
	fname = [titleCase(hdr.patient_name) '_' hdr.patientID '_' findTracer({hdr, 'UnknownTracer'}) '_' titleCase(hdr.study)];
end

i = 1;
while i<=length(fname)
	if (fname(i)>='A' && fname(i)<='Z') || ...
			(fname(i)>='a' && fname(i)<='z') || ...
			(fname(i)>='0' && fname(i)<='9') || ...
			any(fname(i)=='_()@')
		i = i+1;
	else
		fname(i)='';
	end
end
fname = [DCSVars.outdir filesep fname];

if ~DCSVars.overwrite && exist([fname '.mat'],'file')
	i = 1;
	while exist([fname '_' num2str(i) '.mat'],'file')
		i = i+1;
	end
	fname = [fname '_' num2str(i) '.mat'];
else
	fname = [fname '.mat'];
end
		
hdr.DICOMConversionNotes = ['Converted ' datestr(now,'yyyy-mm-dd HH:MM:SS') ' by ' mfilename ' (V' DCSVars.version ')'];
try
	save(fname,'-mat','-v7.3','hdr','vol','infodcm');
	DCSlogmsg(['Saved to:' fname]);
catch
	DCSlogmsg(' - Error while saving');
	savecurrentfname
	return
end
DCSVars.currentfname = ''; % this is no longer the current file to process


%% Switch string to title case
function str = titleCase(str)
str = lower(str);
indx = [0 find(str==' ')];
for i=1:length(indx)
	if length(str)>indx(i)
		str(indx(i)+1) = upper(str(indx(i)+1));
	end
end


%% Checks if user requested to close the server
function status = serverStopRequested(DCSVars)
>>>>>>> f95da97b0f7b59512a505a4b1d7d27ebe583c1e6
status = strcmpi(get(DCSVars.timerh,'Running'),'off');