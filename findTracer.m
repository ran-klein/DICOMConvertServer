<<<<<<< HEAD
% FINDTRACER - determines the tracer from a description string. str may
% be a filename or another describing string. If a tracer cannot be
% determined tracer is returned with an empty value;
%
% tracer = findTracer(strs,delimeters,tracers)
% tracer = findTracer(strs) - determine the tracer from the cell of strings
% strs, whith the first value having highest priority and the last having
% the lovest priority. strs can also contain a valid header structure, and
% then the fields will be scaned in a preset hirarchy.
%
% tracer = findTracer(strs,delimeters) - same as above, but uses a
% specified set of text deilmeters.
%
% tracer = findTracer(strs,delimeters,tracers) - same as above, but uses a
% specified set of tracers
%
% tracer = findTracer(0,delimeters,tracers) - set tracers as the default
% dictionary.
%
% tracer = findTracer(0,delimeters,) - reverts to the default tracer
% dictionary.
%
% See also: findImageType, findState

% Ran Klein sometime in summer 2007
% mmodified:
% 2008-06-01 - RK - addition of support of hdr structure in strs
%                   cells.


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function tracer = findTracer(strs,delimeters,tracers)
	
persistent tracerDictionary
if isempty(ver('MATLAB')) || verLessThan('MATLAB','7.14')
	iptchecknargin(1,3,nargin,mfilename);
else
	narginchk(1,3);
end

if nargin<3
	if isempty(tracerDictionary)
		tracers = {{'Rubidium','Rb'};...
			{'Ammonia','NH3','Npg'};...
			{'FDG','Fluorodeoxyglucose'};...
			{'Fluoride','F'};...
			{'Rolipram'};...
			{'Acetate'};...
			{'HED'};...
			{'Sestamibi','MIBI','Tetrofosmin','Tc','Technetium'};...
			{'Perfusion CT','CT','Perf'};...
			{'Water','OWater','H2O'};...
			{'Flurpiridaz','BMS7470158','BMS74701582'};...
			{'Contrast CT','Omnipaque'}};
	else
		tracers = tracerDictionary;
	end
end

% -------------------------------------------------------------------------
% we are defining the dictionary
if isnumeric(strs) 
	tracerDictionary = tracers;
	tracer = tracerDictionary;
	return
end
% -------------------------------------------------------------------------

tracerlist = cellarray(tracers);
if nargin<2 || isempty(delimeters)
	delimeters = ' ,_\1234567890.^-';
end

if ~iscell(strs)
	strs = {strs};
end

tracer = '';
done = false;
i=1;
while ~done && i<=length(strs)
	if isstruct(strs{i})
		strs = [strs(1:i-1) strs{i}.tracer, strs{i}.study, strs{i}.filename, strs{i}.nativefile, strs(i+1:end)];
	else
		str = strs{i};
		while ~done && ~isempty(str)
			[a, str] = strtok(str,delimeters);
			indx = strmatch(lower(a),lower(tracerlist),'exact');
			if length(indx)==1
				tracer = tracerlist{indx};
				done = true;
			end
		end
		i=i+1;
	end
end

if isempty(tracer)
	tracer = 'Unknown';
else
	done = false;
	i = 1;
	while ~done && i<=length(tracers)
		if any(strcmpi(tracer,tracers{i}))
			done = true;
			tracer = tracers{i}{1};
		else
			i = i+1;
		end
	end
end

%% CELLARRAY - convert a cell of cells into an array of cells
function array = cellarray(cells)
if iscell(cells)
	array = {};
	for i = 1:length(cells)
		array = [array cellarray(cells{i})];
	end
else
	array = cells;
=======
% FINDTRACER - determines the tracer from a description string. str may
% be a filename or another describing string. If a tracer cannot be
% determined tracer is returned with an empty value;
%
% tracer = findTracer(strs,delimeters,tracers)
% tracer = findTracer(strs) - determine the tracer from the cell of strings
% strs, whith the first value having highest priority and the last having
% the lovest priority. strs can also contain a valid header structure, and
% then the fields will be scaned in a preset hirarchy.
%
% tracer = findTracer(strs,delimeters) - same as above, but uses a
% specified set of text deilmeters.
%
% tracer = findTracer(strs,delimeters,tracers) - same as above, but uses a
% specified set of tracers
%
% tracer = findTracer(0,delimeters,tracers) - set tracers as the default
% dictionary.
%
% tracer = findTracer(0,delimeters,) - reverts to the default tracer
% dictionary.
%
% See also: findImageType, findState

% Ran Klein sometime in summer 2007
% mmodified:
% 2008-06-01 - RK - addition of support of hdr structure in strs
%                   cells.


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function tracer = findTracer(strs,delimeters,tracers)
	
persistent tracerDictionary
if isempty(ver('MATLAB')) || verLessThan('MATLAB','7.14')
	iptchecknargin(1,3,nargin,mfilename);
else
	narginchk(1,3);
end

if nargin<3
	if isempty(tracerDictionary)
		tracers = {{'Rubidium','Rb'};...
			{'Ammonia','NH3','Npg'};...
			{'FDG','Fluorodeoxyglucose'};...
			{'Fluoride','F'};...
			{'Rolipram'};...
			{'Acetate'};...
			{'HED'};...
			{'Sestamibi','MIBI','Tetrofosmin','Tc','Technetium'};...
			{'Perfusion CT','CT','Perf'};...
			{'Water','OWater','H2O'};...
			{'Flurpiridaz','BMS7470158','BMS74701582'};...
			{'Contrast CT','Omnipaque'}};
	else
		tracers = tracerDictionary;
	end
end

% -------------------------------------------------------------------------
% we are defining the dictionary
if isnumeric(strs) 
	tracerDictionary = tracers;
	tracer = tracerDictionary;
	return
end
% -------------------------------------------------------------------------

tracerlist = cellarray(tracers);
if nargin<2 || isempty(delimeters)
	delimeters = ' ,_\1234567890.^-';
end

if ~iscell(strs)
	strs = {strs};
end

tracer = '';
done = false;
i=1;
while ~done && i<=length(strs)
	if isstruct(strs{i})
		strs = [strs(1:i-1) strs{i}.tracer, strs{i}.study, strs{i}.filename, strs{i}.nativefile, strs(i+1:end)];
	else
		str = strs{i};
		while ~done && ~isempty(str)
			[a, str] = strtok(str,delimeters);
			indx = strmatch(lower(a),lower(tracerlist),'exact');
			if length(indx)==1
				tracer = tracerlist{indx};
				done = true;
			end
		end
		i=i+1;
	end
end

if isempty(tracer)
	tracer = 'Unknown';
else
	done = false;
	i = 1;
	while ~done && i<=length(tracers)
		if any(strcmpi(tracer,tracers{i}))
			done = true;
			tracer = tracers{i}{1};
		else
			i = i+1;
		end
	end
end

%% CELLARRAY - convert a cell of cells into an array of cells
function array = cellarray(cells)
if iscell(cells)
	array = {};
	for i = 1:length(cells)
		array = [array cellarray(cells{i})];
	end
else
	array = cells;
>>>>>>> f95da97b0f7b59512a505a4b1d7d27ebe583c1e6
end