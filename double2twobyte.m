% DOUBLE2TWOBYTE - converts a dynamic image sequence represented by double
% values to two byte values (1 - 32766) and returns an array to scale each 
% frame.
%
% [img quantFactor] = double2twobyte(img)
%
% Note: This function is used in conjunction with ReorientTool so that
% volumetric image data can be rotated using resectd.dll
% Examples can be found in ProcessFAD.m:
%        [uptake hdr.quant_uptake] = double2twobyte(cMapVol(:,:,:,myofact));
% 		 [opt, vol] = reorientTool(uptake, hdr, opt);
%
% By Ran Klein 16-Oct-2006


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function [img, quantFactor] = double2twobyte(img)

if ndims(img)<=3 % static image
	% 	quantFactor = max(img(:))/32766;
	% 	img = max(1,ceil(img / quantFactor));

	quantFactor = max( max(img(:))/32767,...
		min(img(:))/-32768 );
	if quantFactor == 0
		quantFactor = 1;
	end
	img = round(img / quantFactor);
		
elseif ndims(img)==4 % dynamic sequence
	nfr = size(img,4);
	quantFactor = zeros(nfr,1);
	for i=1:nfr
% 		quantFactor(i) = max(max(max(img(:,:,:,i))))/32766;
% 		img(:,:,:,i) = max(1,ceil(img(:,:,:,i) / quantFactor(i)));
		
		quantFactor(i) = max( max(max(max(img(:,:,:,i))))/32767,...
			min(min(min(img(:,:,:,i))))/-32768 );
		if quantFactor(i)==0
			quantFactor(i) = 1;
		end
		img(:,:,:,i) = round(img(:,:,:,i) / quantFactor(i));
		
	end
	% Equivalent implementation is faster, but uses more memory.
% 	for i=1:nfr
% 		quantFactor(i) = max(max(max(img(:,:,:,i))))/32766;
% 	end
% 	s = size(img);
% 	img = reshape(max(1,ceil(reshape(img,[prod(s(1:3)),s(4)])/quantFactor')) ,s);
end
img = int16(img);