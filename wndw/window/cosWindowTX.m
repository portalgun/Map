function W = cosWindowTX(numPixT, numPixX, dskDmPixT, dskDmPixX, rmpWidthPixT, rmpWidthPixX, bSym,bPLOT)

% function W = cosWindowTX(numPixOrig,dskDmPix,rmpWidthPixX2,bSym,bPLOT)
% 
%   example call: W = cosWindowTX(16, 16, 4, 0,12, 16,1 ,1);
%
%                 W = cosWindowTX(16, 16, 2, 0,12, 16,1 ,1);
%
%                 W = cosWindowTX(16, 16, 0, 0,12, 16,1 ,1);
%
% build cosine window with a flatop
% 
% numPixT:       size of window in T,
% numPixX:       size of window in X,
% dskDmPixT:     t diameter of disk (flattop) in pixels
% dskDmPixX:     x diameter of disk (flattop) in pixels
% rmpWidthPixT:  t diameter of ramp in pixels (i.e. twice the ramp radius)
% rmpWidthPixX:  x diameter of ramp in pixels (i.e. twice the ramp radius)
% bSym:          if numPix is even, boolean to make symmetric flattop cos
%                window profile on either side of zero
% bPLOT:         1 -> plot 
%                0 -> not
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

if dskDmPixT + rmpWidthPixT > numPixT
   disp(['cosWindowTX: WARNING! disk + ramp radius exceeds image size along T']);
end
if dskDmPixX + rmpWidthPixX > numPixX
   disp(['cosWindowTX: WARNING! disk + ramp radius exceeds image size along X']);
end

if ~exist('bSym','var') || isempty(bSym)
    bSym = 0;
end
if ~exist('bPLOT','var') || isempty(bPLOT)
   bPLOT = 0; 
end


% CONVERT DIAMETER TO RADIUS
dskRadiusPixT = dskDmPixT/2;
dskRadiusPixX = dskDmPixX/2;
rmpRadiusPixT = rmpWidthPixT/2;
rmpRadiusPixX = rmpWidthPixX/2;

% % % MESHGRID LOCATIONS
% % if bSym == 0
% %     [T X] = meshgrid(Wave.smpPos(1,numPix));
% % elseif bSym == 1 && mod(numPix,2)==0
% %     x = Wave.smpPos(1,numPix);
% %     x = x+diff(x(1:2))/2;
% %     [T X] = meshgrid(x); 
% % else
% %     error(['cosWindowTX: WARNING! unhandled bSym and numPix values']);
% % end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD FLATTOP COSINE WINDOW in X %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RX =  Wave.smpPos(1,numPixX);

X = ones(numPixX,1);
% MAKE RAMP
freqcpp = 1./(2*rmpRadiusPixX); % cycles per pixel
X(abs(RX)>dskRadiusPixX) = 0.5.*(1 + cos(2.*pi.*freqcpp*(abs(RX(abs(RX)>dskRadiusPixX))-dskRadiusPixX)));

% SET VALUES OUTSIDE OF RAMP TO ZERO
X(abs(RX)>(dskRadiusPixX+rmpRadiusPixX))= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD FLATTOP COSINE WINDOW in T %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RT =  Wave.smpPos(1,numPixT);

T = ones(numPixT,1);
% MAKE RAMP 
freqcpp = 1./(2*rmpRadiusPixT); % cycles per pixel
T(abs(RT)>dskRadiusPixT) = 0.5.*(1 + cos(2.*pi.*freqcpp*(abs(RT(abs(RT)>dskRadiusPixT))-dskRadiusPixT)));

% SET VALUES OUTSIDE OF RAMP TO ZERO
T(abs(RT)>(dskRadiusPixT+rmpRadiusPixT))= 0;


W = T*X';
if bPLOT
   figure('position',[680   666   805   368]); 
   subplot(1,3,1);
   imagesc(RX,RT,W);
   axis square
   axis xy
   Fig.format(['X'],['T'],['DskT=' num2str(dskDmPixT) '; DskX=' num2str(dskDmPixX)]);
   
   ind = floor(size(W,1)./2 + 1);
   subplot(1,3,2);
   plot(RX,W(ind,:),'k');
   axis square
   axis xy
   Fig.format(['X'],['W'],['RmpPixX=' num2str(rmpWidthPixX)]);
   
   subplot(1,3,3);
   plot(RT,W(:,ind),'k');
   axis square
   axis xy
   Fig.format(['X'],['T'],['RmpPixT=' num2str(rmpWidthPixT)]);
   
end
