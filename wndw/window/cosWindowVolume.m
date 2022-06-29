function W = cosWindowVolume(WszRCT,Wk,bPLOT)

% function W = cosWindowVolume(WszRCT,Wk,bPLOT)
% 
%   example calls: W = cosWindowVolume([26 26 26],1,1);
%
% build volume cosine window with a radius equal to the smallest of the
% dimensions. cos windows should typically be used only with square tensors
%  
% WszRCT:     size of patch (row,col,frame) in pixels    [ 1 x 3 ]
%            if scalar entered, function converts to     [ 1 x 3 ]
% Wk:        amount of the patch to cosine window
%            1, .5, .25, .125, .0625  are the only acceptable values
% bPLOT:     1 -> plot, 
%            0 -> don't (default == 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W:         cosine window                              [ R x C x F ]

if ~exist('Wk','var') || isempty(Wk) Wk    = 1; end
if ~exist('bPLOT','var')             bPLOT = 0; end
if length(WszRCT) == 1   WszRCT = repmat(WszRCT,1,2); end
if length(WszRCT) == 3 && (WszRCT(1) == 1 || WszRCT(2) == 1)
    bSlice1D = 1;
    WszRCT = max(WszRCT)*[1 1 1];
else
    bSlice1D = 0;
end
WszMin = min(WszRCT);

% CHECK FOR LEGAL VALUES OF Wk (poorly written)
if rem(1,Wk) ~= 0
   error(['cosWindowVolume: Wk must divide zero evenly. Current value = ' num2str(Wk)]);
end

% POSITION SAMPLES
if     mod(WszRCT(2),2) == 0 ||  WszRCT(2)==1, Xdeg = Wave.smpPos(WszRCT(2),WszRCT(2));
elseif mod(WszRCT(2),2) == 1,                  Xdeg = [(-(WszRCT(2)-1)/2):((WszRCT(2)-1)/2)]./(WszRCT(2)-1);
else   error(['cosWindowVolume: WARNING! WszRCT must be integer valued!!!']);
end
if     mod(WszRCT(1),2) == 0 ||  WszRCT(1)==1, Ydeg = Wave.smpPos(WszRCT(1),WszRCT(1));
elseif mod(WszRCT(1),2) == 1,                  Ydeg = [(-(WszRCT(1)-1)/2):((WszRCT(1)-1)/2)]./(WszRCT(1)-1);
else   error(['cosWindowVolume: WARNING! WszRCT must be integer valued!!!']);
end
if     mod(WszRCT(3),2) == 0 ||  WszRCT(3)==1, Zdeg = Wave.smpPos(WszRCT(3),WszRCT(3));
elseif mod(WszRCT(3),2) == 1,                  Zdeg = [(-(WszRCT(3)-1)/2):((WszRCT(3)-1)/2)]./(WszRCT(3)-1);
else   error(['cosWindowVolume: WARNING! WszRCT must be integer valued!!!']);
end

% BUILD COSINE 
[W,Xdeg,Ydeg,Zdeg] = cosdRadialVolume(Xdeg,Ydeg,Zdeg,1/Wk,.5,.5,0);

% RADIAL DISTANCE FROM CENTER
Rdeg = sqrt(Xdeg.^2 + Ydeg.^2 + Zdeg.^2);
% ZERO RESULTS
W(Rdeg>=Wk./2) = 0;

if bSlice1D == 1
   W = W(floor(1+size(W,1)/2),:,:);
end

if bPLOT
    if bSlice1D == 0
        figure;
        for i = 1:size(Zdeg,3)
            imagesc(Xdeg(1,:,i),Ydeg(:,1,i)',W(:,:,i));
            Fig.format('X (deg)','Y (deg)',['Z=' num2str(Zdeg(1,1,i),'%.2f')]);
            minmax(W);
            caxis([minmax(W)]);
            axis square;
            % colormap gray
            pause(.2);
        end
    else
        figure;
        for i = 1:size(Zdeg,3)
        plot(unique(Xdeg),squeeze(W(floor(size(W,1)/2+1),:,i)),'k-','linewidth',2);
        Fig.format('X (deg)','W',['Z=' num2str(Zdeg(1,1,i),'%.2f')]);
        axis square
        ylim([minmax(W)]);
        pause(.2); 
        end
    end
end