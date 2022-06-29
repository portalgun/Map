function W = cosWindowXYT(PszXYT, rmpDmXYT, dskDmXYT, bSymS,  bPLOT)

% function W = cosWindowXYT(PszXYT, rmpDmXYT, dskDmXYT, bSymS, bPLOT)
% 
%   example call: W = cosWindowXYT([16 16 16], [16 16 8], [0 0 8], 1, 1);
%
% build x,y,t cosine window with a flatop
% 
% PszXYT:       x,y,t size of window in pixels                                [ 1 x 3 ] 
% dskDmXYT:     x,y,t diameter of disk (flattop) in pixels                    [ 1 x 3 ] 
% rmpDmXYT:     x,y,t diameter of ramp in pixels (i.e. twice the ramp radius) [ 1 x 3 ] 
% bSymS:        boolean indicating whether space components should be 
%               symmetric or not. if not, WXY computed via outer product 
% bPLOT:        1 -> plot 
%               0 -> not
% %%%%%%%%%%%%%%%%%%%
% W:           x,y,t cosine window with arbitrary disk and ramp width
%       ______
%      /      \
%     /        \
% ___/          \___
%


% INPUT CHECKING
if ~exist('rmpDmXYT','var') || isempty(rmpDmXYT) rmpDmXYT = PszXYT;          end
if ~exist('dskDmXYT','var') || isempty(dskDmXYT) dskDmXYT = PszXYT-rmpDmXYT; end
if ~exist('bSymS','var')    || isempty(bSymS)    bSymS    = 1;               end
if ~exist('bPLOT','var')    || isempty(bPLOT),   bPLOT    = 0;               end

% INPUT HANDLING
if bSymS == 1 && (dskDmXYT(1) ~= dskDmXYT(2) || rmpDmXYT(1) ~= rmpDmXYT(2) )
   error(['cosWindowXYT: WARNING! parameters inconsistent with spatial symmmetry. rmpDmXYT=[' num2str(rmpDmXYT) '] and dskDmXYT=[' num2str(dskDmXYT) ']']);
end

numPixX      =   PszXYT(1); numPixY      =   PszXYT(2);      numPixT =   PszXYT(3);
dskDmPixX    = dskDmXYT(1); dskDmPixY    = dskDmXYT(2);    dskDmPixT = dskDmXYT(3);
rmpWidthPixX = rmpDmXYT(1); rmpWidthPixY = rmpDmXYT(2); rmpWidthPixT = rmpDmXYT(3);
if dskDmPixX + rmpWidthPixX > numPixX, disp(['cosWindowXYT: WARNING! disk + ramp widths exceeds image size along WX']); end
if dskDmPixY + rmpWidthPixY > numPixY, disp(['cosWindowXYT: WARNING! disk + ramp widths exceeds image size along WY']); end
if dskDmPixT + rmpWidthPixT > numPixT, disp(['cosWindowXYT: WARNING! disk + ramp widths exceeds image size along WT']); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD SYMMETRIC FLATTOP COSINE WINDOW in WXY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bSymS == 1 

[X Y]= meshgrid(Wave.smpPos(1,numPixX),Wave.smpPos(1,numPixY));

if mod(numPixX,2) == 0
   X = X+diff(X(1,1:2))/2;  
end

if mod(numPixY,2) == 0
   Y = Y+diff(Y(1:2,1))/2;  
end

RXY =  sqrt(X.^2+Y.^2);

% CONVERT DIAMETER TO RADIUS
dskRadiusPix =    dskDmPixX/2; 
rmpRadiusPix = rmpWidthPixX/2; 
% MAKE RAMP
freqcpp = 1./(2*rmpRadiusPix); % cycles per pixel
WXY = ones(size(RXY));
WXY(RXY>dskRadiusPix) = 0.5.*(1 + cos(2.*pi.*freqcpp*(RXY(RXY>dskRadiusPix)-dskRadiusPix)));
% SET VALUES OUTSIDE OF RAMP TO ZERO
WXY(RXY>(dskRadiusPix+rmpRadiusPix)) = 0;

elseif (bSymS == 0)
[X Y]=meshgrid(Wave.smpPos(1,numPixX),Wave.smpPos(1,numPixY));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD FLATTOP COSINE WINDOW in WX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RX =  Wave.smpPos(1,numPixX);
if mod(numel(RX),2) == 0
   RX = RX+diff(RX(1:2))/2;  
end
RX;
% CONVERT DIAMETER TO RADIUS
dskRadiusPixX =    dskDmPixX/2; 
rmpRadiusPixX = rmpWidthPixX/2; 
% MAKE RAMP
freqcpp = 1./(2*rmpRadiusPixX); % cycles per pixel
WX = ones(numPixX,1);
WX(abs(RX)>dskRadiusPixX) = 0.5.*(1 + cos(2.*pi.*freqcpp*(abs(RX(abs(RX)>dskRadiusPixX))-dskRadiusPixX)));
% SET VALUES OUTSIDE OF RAMP TO ZERO
WX(abs(RX)>(dskRadiusPixX+rmpRadiusPixX)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD FLATTOP COSINE WINDOW in WY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RY =  Wave.smpPos(1,numPixY);
if mod(numel(RY),2) == 0
   RY = RY+diff(RY(1:2))/2;  
end
RY;
% CONVERT DIAMETER TO RADIUS
dskRadiusPixY =    dskDmPixY/2; 
rmpRadiusPixY = rmpWidthPixY/2; 
% MAKE RAMP
freqcpp = 1./(2*rmpRadiusPixY); % cycles per pixel
WY = ones(numPixY,1);
WY(abs(RY)>dskRadiusPixY) = 0.5.*(1 + cos(2.*pi.*freqcpp*(abs(RY(abs(RY)>dskRadiusPixY))-dskRadiusPixY)));
% SET VALUES OUTSIDE OF RAMP TO ZERO
WY(abs(RY)>(dskRadiusPixY+rmpRadiusPixY)) = 0;

WXY = WY*WX';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD FLATTOP COSINE WINDOW in WT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RT =  Wave.smpPos(1,numPixT);
if mod(numel(RT),2) == 0
   RT = RT+diff(RT(1:2))/2;
end
% CONVERT DIAMETER TO RADIUS
dskRadiusPixT =    dskDmPixT/2;
rmpRadiusPixT = rmpWidthPixT/2;
% MAKE RAMP 
freqcpp = 1./(2*rmpRadiusPixT); % cycles per pixel
WT = ones(numPixT,1);
WT(abs(RT)>dskRadiusPixT) = 0.5.*(1 + cos(2.*pi.*freqcpp*(abs(RT(abs(RT)>dskRadiusPixT))-dskRadiusPixT)));
% SET VALUES OUTSIDE OF RAMP TO ZERO
WT(abs(RT)>(dskRadiusPixT+rmpRadiusPixT))= 0;

W = bsxfun(@times,WXY,permute(WT,[3 2 1]));
if bPLOT
   figure('position',[680   666   805   368]); 
   for i = 1:size(W,3)
       clf;
       subplot(1,3,1);
       imagesc(X(1,:),Y(:,1)',W(:,:,i));
       axis square
       axis xy
       Fig.format(['WX'],['WY'],['RmpT=' num2str(rmpDmXYT(3)) ', DskT=' num2str(dskDmPixT)]);
       caxis([0 1])

       ind = floor(size(W,1)./2 + 1);
       subplot(1,3,2);
       plot(X(1,:),W(ind,:,i),'k');
       axis square
       axis xy
       Fig.format(['WX'],['WY'],['RmpPixX=' num2str(rmpWidthPixX)]);
       ylim([0 1])

       subplot(1,3,3);
       plot(Y(:,1)',W(:,ind,i),'k');
       axis square
       axis xy
       Fig.format(['WY'],['W'],['RmpPixY=' num2str(rmpWidthPixY)]);
       ylim([0 1])

       sgtitle(['t=' num2str(i)])
       pause(.05); 
   end
end
