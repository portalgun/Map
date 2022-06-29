function WspaceTime = cosWindowTemporal(numSmpSpaceXY,numSmpTime,numSensors,winFactor)

% function WspaceTime = cosWindowTemporal(numSmpSpaceXY,numSmpTime,winFactor)
% 
%   example calls: WspaceTime = cosWindowTemporal([128 128],16,1);
%                  WspaceTime = cosWindowTemporal([128 128],16,1,.5);
%
% build 3D cosine window with a radius that is equal to the smaller of the two
% dimensions. it is recommended that cos windows are used only with square
% patches.
%  
% numSmpSpaceXY:   [1x2] vector containing the size of the patch in pixels
% numSmpTime:      scalar containing the number of temporal samples
% numSensors:      number of spatial sensor types
%                  1 -> luminance sensors only, for example
%                  3 -> rgb       sensors
% winFactor:       amount of the patch to cosine window
%                  1, .5, .25, .125, .0625  are the only exceptable values
% bPLOT:           1 -> plot, 0 -> don't (default == 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cosWin:  cosine window

if ~exist('winFactor','var') || isempty(winFactor)
    winFactor = [];
end
% FULL COSINE WINDOW (SPACE)
Wspace              = cosWindow(numSmpSpaceXY,winFactor);
% HALF-COSINE WINDOW (TIME)
Wtime               = cosWindow(numSmpTime*2*[1 1],winFactor);
Wtime               = Wtime(numSmpTime+1,:);
Wtime(1:numSmpTime) = [];
% SPACE-TIME WINDOW
WspaceTime          = bsxfun(@times,Wspace,permute(Wtime,[4 3 2 1]));
WspaceTime          = reshape(WspaceTime,[numSmpSpaceXY numSensors numSmpTime]);