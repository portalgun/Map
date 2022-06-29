function W = cosWindow(WszRCT,Wk,bPLOT)

% function W = cosWindow(WszRCT,Wk,bPLOT)
% 
%   example calls: W = cosWindow(size(I));
%                  W = cosWindow(size(I),.5);
%
% build 2D cosine window with a radius that is equal to the smaller of the two
% dimensions. it is recommended that cos windows are used only with square
% patches.
%  
% WszRCT:    vector indicating row,col,time size of patch in pixels 
%            if size(WszRCT) = [ 1 x 1 ], function converts to   [ 1 x 2 ]
%            if size(WszRCT) = [ 1 x 2 ] or [ 2 x 1 ] AND
%            input size is WszRCT = [ 1  N ], creates a 1D slice
% Wk:        amount of the patch to cosine window
%            1, .5, .25, .125, .0625  are the only acceptable values
% bPLOT:     1 -> plot, 
%            0 -> don't (default == 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W:  cosine window

if ~exist('Wk','var') || isempty(Wk) Wk    = 1; end
if ~exist('bPLOT','var')             bPLOT = 0; end
if length(WszRCT) == 1 WszRCT = repmat(WszRCT,1,2); end
if (length(WszRCT) == 2)  && min(WszRCT) == 1
    bSlice1D = 1;
    WszRCT = max(WszRCT)*[1 1];
else
    bSlice1D = 0;
end
if length(WszRCT) == 3
   W = cosWindowVolume(WszRCT,Wk,bPLOT); 
   return;
end
WszMin = min(WszRCT);

% CHECK FOR LEGAL VALUES OF Wk (poorly written)
if rem(1,Wk) ~= 0
   error(['cosWindow: Wk must divide zero evenly. Current value = ' num2str(Wk)]);
end

% POSITION SAMPLES
if mod(WszMin,2) == 0
    pos = Wave.smpPos(WszMin,WszMin); % [(-WszMin/2):(WszMin/2-1)]./WszMin;
elseif mod(WszMin,2) == 1
    % NOTE! DO NOT USE Wave.smpPos.m HERE
    pos = [(-(WszMin-1)/2):((WszMin-1)/2)]./(WszMin-1);
else
    error('cosWindow: WHOA NELLY!');
end
[xx,yy] = meshgrid(pos);
rr = sqrt(xx.^2 + yy.^2);
W = 0.5.*(1 + cos(2.*pi.*rr./Wk));
W(rr>=Wk./2) = 0;

if bSlice1D == 1
   W = W(floor(1+size(W,1)/2),:);
end

if bPLOT
    if bSlice1D == 0
    figure; imagesc(W); 
    axis tight; box on; zlim([0 1]); axis image;
    else
    figure; 
    plot(unique(xx),W(floor(size(W,1)/2+1),:),'k-','linewidth',2);
    end
end