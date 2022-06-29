function W = cosWindowXT(PszXT, dskDmXT, rmpDmXT, bSym,bPLOT)

% function W = cosWindowXT(PszXT, dskDmXT, rmpDmXT, bSym,bPLOT)
% 
%   example call: W = cosWindowXT([16 16], [0 4], [16 12], 1 ,1);
%
%                 W = cosWindowXT([16 16], [0 2], [16 14], 1 ,1);
%
%                 W = cosWindowXT([16 16], [0 0], [16 16], 1 ,1);
%
% build cosine window with a flatop
% 
% PszXT:       x,t size of window in pixels                                [ 1 x 2 ] 
% dskDmXT:     x,t diameter of disk (flattop) in pixels                    [ 1 x 2 ] 
% rmpDmXT:     x,t diameter of ramp in pixels (i.e. twice the ramp radius) [ 1 x 2 ] 
% bSym:        if numPix is even, boolean to make symmetric flattop cos
%              window profile on either side of zero
%              1 -> yes
%              0 -> no
% bPLOT:       1 -> plot 
%              0 -> not
% %%%%%%%%%%%%%%%%%%%
% W:             x,t cosine window with arbitrary disk and ramp width
%       ______
%      /      \
%     /        \
% ___/          \___
%

% if length(numPix)==1
%    numPix = numPix*ones(1,2); 
% end


numPixX = PszXT(1);
numPixT = PszXT(2);

dskDmPixX = dskDmXT(1);
dskDmPixT = dskDmXT(2);

rmpWidthPixX= rmpDmXT(1);
rmpWidthPixT= rmpDmXT(2);

if dskDmPixT + rmpWidthPixT > numPixT
   disp(['cosWindowXT: WARNING! disk + ramp radius exceeds image size along WT']);
end
if dskDmPixX + rmpWidthPixX > numPixX
   disp(['cosWindowXT: WARNING! disk + ramp radius exceeds image size along WX']);
end

if ~exist('bSym','var') || isempty(bSym)
    bSym = 0;
end
if ~exist('bPLOT','var') || isempty(bPLOT)
   bPLOT = 0; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD FLATTOP COSINE WINDOW in WX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RX =  Wave.smpPos(1,numPixX);
if bSym == 1
   RX = RX+diff(RX(1:2))/2;  
end
RX;
% CONVERT DIAMETER TO RADIUS
dskRadiusPixX = dskDmPixX/2;
dskRadiusPixT = dskDmPixT/2;

rmpRadiusPixX = rmpWidthPixX/2;
rmpRadiusPixT = rmpWidthPixT/2;


WX = ones(numPixX,1);
% MAKE RAMP
freqcpp = 1./(2*rmpRadiusPixX); % cycles per pixel
WX(abs(RX)>dskRadiusPixX) = 0.5.*(1 + cos(2.*pi.*freqcpp*(abs(RX(abs(RX)>dskRadiusPixX))-dskRadiusPixX)));
% SET VALUES OUTSIDE OF RAMP TO ZERO
WX(abs(RX)>(dskRadiusPixX+rmpRadiusPixX)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD FLATTOP COSINE WINDOW in WT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RT =  Wave.smpPos(1,numPixT);
if bSym == 1
   RT = RT+diff(RT(1:2))/2;
end
WT = ones(numPixT,1);
% MAKE RAMP 
freqcpp = 1./(2*rmpRadiusPixT); % cycles per pixel
WT(abs(RT)>dskRadiusPixT) = 0.5.*(1 + cos(2.*pi.*freqcpp*(abs(RT(abs(RT)>dskRadiusPixT))-dskRadiusPixT)));
% SET VALUES OUTSIDE OF RAMP TO ZERO
WT(abs(RT)>(dskRadiusPixT+rmpRadiusPixT))= 0;

W = WT*WX';
if bPLOT
   figure('position',[680   666   805   368]); 
   subplot(1,3,1);
   imagesc(RX,RT,W);
   axis square
   axis xy
   Fig.format(['WX'],['WT'],['DskT=' num2str(dskDmPixT) '; DskX=' num2str(dskDmPixX)]);
   
   ind = floor(size(W,1)./2 + 1);
   subplot(1,3,2);
   plot(RX,W(ind,:),'k');
   axis square
   axis xy
   Fig.format(['WX'],['W'],['RmpPixX=' num2str(rmpWidthPixX)]);
   
   subplot(1,3,3);
   plot(RT,W(:,ind),'k');
   axis square
   axis xy
   Fig.format(['WT'],['W'],['RmpPixT=' num2str(rmpWidthPixT)]);
   
end
