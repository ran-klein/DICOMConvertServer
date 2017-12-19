% DCSlogmsg - Log a message to the log GUI and to log file.
%
% Usage: fig1 = DCSlogmsg(str,flag)
% Input parameters:
% -----------------
% str - the string message
% flag - message display format:
%          2 - replace current line
% 		   1 - append to end of current line
% 	    else - new line
%
% Output parameters
% -----------------
% fig1 - figure handle of the output GUI.
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


function fig1 = DCSlogmsg(str,flag)
global DCSVars
persistent fig
if nargin<2
	flag = 0;
end
fig1 = fig;

%% Create the figure
if isempty(fig) || ~ishandle(fig)
	fig = findobj(get(0,'children'),'flat','name','DICOM Convert Server'); % recycle window
	if isempty(fig)
		fig = figure;
	end
	set(fig,'name','DICOM Convert Server','numbertitle','off',...
		'CloseRequestFcn', 'dicomConvertServer(''stop'');',... 
		'DeleteFcn',@TerminateLog,...
		'Units','normalized','Position',[0.2 0.2 0.6 0.6],'menubar','none','UserData','Running');
	movegui(fig,'center');
	texth = uicontrol(fig,'tag','TextBox','style','listbox','Units','normalized','Position',[0.05 0.05 0.9 0.9],...
		'enable','inactive','horizontalalignment','left',...
		'fontname','Courier New','fontsize',8);
	setappdata(fig,'TextBoxHandle',texth);
	strs = {'Starting DICOM Convert Server Log'};
else
	texth = getappdata(fig,'TextBoxHandle');
	strs = get(texth,'string');
end

%% Add a log message
if nargin>=1
	if flag ==2 % replace line
		strs{end} = [datestr(clock) ' : ' str];
	elseif flag == 1 % append to line
		strs{end} = [strs{end},str];
	else % new line
		strs{length(strs)+1} = [datestr(clock) ' : ' str];
	end

	set(texth,'units','points');
	pos = get(texth,'position');
	set(texth,'units','normalized');
	while 1.3*length(strs)*get(texth,'fontsize')>pos(4)
		strs = strs(2:end);
	end
	set(texth,'String',strs);
	drawnow;

	try
		if flag==0
			fprintf(DCSVars.logfileh,'\n%s',[datestr(clock) ' : ' str]);
		elseif flag==1
			fprintf(DCSVars.logfileh,'%s',str);
		end
	catch
	end
end



function TerminateLog(hObject, eventdata)
global DCSVars
try	fclose(DCSVars.logfileh);	catch, end
