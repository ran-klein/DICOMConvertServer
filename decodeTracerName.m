% decodeTraceName - Looks up the tracer name using code scheme designators.
%
% Usage:
% tracer = decodeTracerName(dcmInfo) - where decodeTracerName is the dicom
% header sub structure with the tracer information.
%
% See also: Radiopharmaceuticals.mat, hdrinitdcm

% By Ran Klein 2013-02-05


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function tracer = decodeTracerName(dcmInfo)

tracer = 'Unknown';

if exist('Radiopharmaceuticals.mat','file')
	load('Radiopharmaceuticals.mat')
	
	if isfield(dcmInfo,'CodeValue')
		i = find(strcmpi(tracerList(:,3),dcmInfo.CodeValue));
		if ~isempty(i)
			tracer = tracerList{i,5};
			if isempty(tracer)
				tracer = tracerList{i,6};
			end
			if isempty(tracer)
				tracer = tracerList{i,7};
			end
		else
			i = find(strcmpi(isotopeList(:,3),dcmInfo.CodeValue));
			if ~isempty(i)
				tracer = isotopeList{i,5};
				if isempty(tracer)
					tracer = tracerList{i,4};
				end
			end
		end
	end
end

if strcmpi(tracer,'Unknown')
	if isfield(dcmInfo,'CodeMeaning')
		tracer = dcmInfo.CodeMeaning;
	else
		switch dcmInfo.CodingSchemeDesignator
			case {'SNM3','99SDM'} % CID 4020 - PET RadioNuclide
				switch dcmInfo.CodeValue
					case 'C-111A1', tracer = '^18^Fluorine';
					case 'C-B1031', tracer = 'Fluorodeoxyglucose F^18^';
					case 'C-B1032', tracer = 'Sodium fluoride F^18^';
					case 'C-159A2', tracer = 'Rb^82^ [^82^Rubidium]';
					case 'C-107A1', tracer = 'N^13^ [^13^Nitrogen]';
					case 'C-105A1', tracer = 'C^11^ [^11^Carbon]';
					case 'C-128A2', tracer = 'Ge^68^ [^68^Germanium]';
					case 'C-155A1', tracer = 'Na^22^ [^22^Sodium]';
				end
			case 'SRT' % CID 4021 - PET Radiopharmacuticals
				switch dcmInfo.CodeValue
					case 'C-B1043', tracer = 'Acetate C^11^';
					case 'C-B103C', tracer = 'Ammonia N^13^';
					case 'C-B103B', tracer = 'Carbon dioxide O^15^';
					case 'C-B1045', tracer = 'Carbon monoxide C^11^';
					case 'C-B103A', tracer = 'Carbon monoxide O^15^';
					case 'C-B103F', tracer = 'Carfentanil C^11^';
					case 'C-B1034', tracer = 'Fluoro-L-dopa F^18^';
					case 'C-B1046', tracer = 'Germanium Ge^68^';
					case 'C-B103D', tracer = 'Glutamate N^13^';
					case 'C-B103E', tracer = 'Methionine C^11^';
					case 'C-B1038', tracer = 'Oxygen O^15^';
					case 'C-B1039', tracer = 'Oxygen-water O^15^';
					case 'C-B1044', tracer = 'Palmitate C^11^';
					case 'C-B1042', tracer = 'Raclopride C^11^';
					case 'C-B1037', tracer = 'Rubidium chloride Rb^82^';
					case 'C-B1047', tracer = 'Sodium Na^22^';
					case 'C-B1033', tracer = 'Spiperone F^18^';
					case 'C-B1036', tracer = 'Thymidine (FLT)F^18^';
				end
		end
	end
end