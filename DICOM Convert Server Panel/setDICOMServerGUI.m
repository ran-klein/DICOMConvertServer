% setDICOMServerGUI - launch a GUI for cnfiguring ans starting the DICOM
% Convert Server.
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

function varargout = setDICOMServerGUI(varargin)
% SETDICOMSERVERGUI M-file for setDICOMServerGUI.fig
%      SETDICOMSERVERGUI, by itself, creates a new SETDICOMSERVERGUI or raises the existing
%      singleton*.
%
%      H = SETDICOMSERVERGUI returns the handle to a new SETDICOMSERVERGUI or the handle to
%      the existing singleton*.
%
%      SETDICOMSERVERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETDICOMSERVERGUI.M with the given input arguments.
%
%      SETDICOMSERVERGUI('Property','Value',...) creates a new SETDICOMSERVERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before setDICOMServerGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to setDICOMServerGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setDICOMServerGUI

% Last Modified by GUIDE v2.5 23-May-2008 23:35:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @setDICOMServerGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @setDICOMServerGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

% Are parameters passed to FlowQuant?
params = nargin && ischar(varargin{1}) && (any(varargin{1}=='\') || any(varargin{1}=='/') || any(varargin{1}(2)==':'));
			   
if nargin && ischar(varargin{1}) && ~params
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before setDICOMServerGUI is made visible.
function setDICOMServerGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

if ~isempty(varargin)
	if iscell(varargin)
		fname = varargin{1};
	else
		fname = varargin;
	end
	if exist(fname,'file')==7 % fname is a directory
		if fname(end) ~= filesep
			fname = [fname filesep];
		end
		fname = [fname 'DICOMServer.dat'];
	end
	setappdata(handles.figure,'SettingsFilename',fname);
end

% Default settings
indir = [rootdir filesep 'DICOMImport'];
outdir = [rootdir filesep 'DataStash'];
logdir = [outdir filesep 'Log'];
tempdir = [outdir filesep 'temp'];
cleanup = true;
overwrite = false;

if exist(filename(handles),'file')
	f = fopen(filename(handles),'r');
	a = '';
	while ~isequal(a,-1)
		a = fgetl(f);
		if ~isequal(a,-1) && ~isempty(a) && a(1)~='%' % not a comment or empty line
			i = find(a=='=');
			if length(i)==1 % valid lines only have one
				temp = strtrim(a(i+1:end));
				if ~isempty(temp) % value is not empty
					if ~isempty(str2num(temp)) % value is a number
						eval([strtrim(a(1:i-1)) ' = ' temp ';']);
					else % value is a string
						eval([strtrim(a(1:i-1)) ' = ''' temp ''';']);
					end
				end % empty value
			end
		end % comment check
	end
	fclose(f);
end

set(handles.IndirEdit,'string',indir);
set(handles.TempdirEdit,'string',tempdir);
set(handles.OutdirEdit,'string',outdir);
set(handles.LogdirEdit,'string',logdir);
% set(handles.CleanupCheckbox,'value',cleanup);
set(handles.OverwriteCheckbox,'value',overwrite);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes setDICOMServerGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = setDICOMServerGUI_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function IndirEdit_Callback(hObject, eventdata, handles)
indir = get(handles.IndirEdit,'string');
if length(indir)<=3
	errordlg({['The selected incoming directory (' indir ') is not valid.'],'Select a new directory.'},'Invalid Incoming Directory','modal');
	set(handles.IndirEdit,'string',[rootdir filesep 'DICOMImport']);
else
	beep
	warndlg({['All files in the incoming directory (' indir ') will be permenently DELETED as part of the conversion.'],'','Please ensure that this is the correct directory.'},'Incoming Directory','modal');
end

% --- Executes during object creation, after setting all properties.
function IndirEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OutdirEdit_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.

function OutdirEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OKButton.
function OKButton_Callback(hObject, eventdata, handles)
saveSettings(handles);

close(handles.figure);


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
close(handles.figure);

% --- Executes on button press in BrowseIndir.
function BrowseIndir_Callback(hObject, eventdata, handles)
path = uigetdir(get(handles.IndirEdit,'string'),'Input directory');
if ~isequal(path,0)
	set(handles.IndirEdit,'String',path);
end
IndirEdit_Callback(hObject, eventdata, handles);


% --- Executes on button press in BrowseOutdir.
function BrowseOutdir_Callback(hObject, eventdata, handles)
path = uigetdir(get(handles.OutdirEdit,'string'),'Output directory');
if ~isequal(path,0)
	set(handles.OutdirEdit,'String',path);
end



function TempdirEdit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function TempdirEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BrowseTempdir.
function BrowseTempdir_Callback(hObject, eventdata, handles)
path = uigetdir(get(handles.TempdirEdit,'string'),'Temporary directory');
if ~isequal(path,0)
	set(handles.TempdirEdit,'String',path);
end


function LogdirEdit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function LogdirEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BrowseLogdir.
function BrowseLogdir_Callback(hObject, eventdata, handles)
path = uigetdir(get(handles.LogdirEdit,'string'),'Logging directory');
if ~isequal(path,0)
	set(handles.LogdirEdit,'String',path);
end



function fname = filename(handles)
fname = getappdata(handles.figure,'SettingsFilename');
if isempty(fname)
	fname = [fileparts(mfilename('fullpath')) filesep 'DICOMServer.dat'];
end



% --- Executes on button press in StartServerButton.
function StartServerButton_Callback(hObject, eventdata, handles)
saveSettings(handles);
dicomConvertServer('start',filename(handles));


function saveSettings(handles)
f = fopen(filename(handles),'w');
fprintf(f,['%% Written by ' strrep(mfilename,'\','\\') '\n' ...
	'%% Date: ' datestr(now) '\n',...
	'indir = ' strrep(get(handles.IndirEdit,'String'),'\','\\') '\n',...
	'outdir = ' strrep(get(handles.OutdirEdit,'String'),'\','\\') '\n',...
	'tempdir = ' strrep(get(handles.TempdirEdit,'String'),'\','\\') '\n',...
	'logdir = ' strrep(get(handles.LogdirEdit,'String'),'\','\\') '\n',...
	'cleanup = 1\n',...
	'overwrite = ' num2str(get(handles.OverwriteCheckbox,'value'))]);
fclose(f);


function OverwriteCheckbox_Callback(hObject, eventdata, handles)
