function W = cosWindowFlattopSpaceTime(PszRC,dskDmSpc,rmpWidthSmpX2spc,bSymSpc,PszT,dskDmTime,rmpWidthSmpX2time,bSymTime,bPLOT)

% function W = cosWindowFlattopSpaceTime(PszRC,dskDmSpc,rmpWidthSmpX2spc,bSymSpc,PszT,dskDmTime,rmpWidthSmpX2time,bSymTime,bPLOT)
% 
%   example call: W = cosWindowFlattopSpaceTime([128 128],0,128,0,16,10,6,1);
%                 
% build cosine window with a flattop in space and time
% 
% PszRC:          size of window in X, Y
% dskDmSpc:      diameter of disk (flattop) in pixels
% rmpWidthSmpX2: diameter of ramp in pixels (i.e. twice the ramp radius)
% bSymSpc:       if numPix is even, boolean to make symmetric flattop cos
%                window profile on either side of zero
% 
%                ADD INPUT PARAM DESCRIPTION
% %%%%%%%%%%%%%%%%%%%
% W:             flattop cosine window
%       ______
%      /      \
%     /        \
% ___/          \___
%

% if length(numPix)==1
%    numPix = numPix*ones(1,2); 
% end

numPix = max(PszRC); 

if dskDmSpc + rmpWidthSmpX2spc > numPix(1)
   disp(['cosWindowFlattopSpaceTime: WARNING! disk + ramp radius exceeds image size']);
end
if ~exist('bPLOT','var') || isempty(bPLOT)
   bPLOT = 0; 
end

Wspc  = cosWindowFlattop(PszRC, dskDmSpc, rmpWidthSmpX2spc, bSymSpc);
Wtime = cosWindowFlattop(PszT ,dskDmTime,rmpWidthSmpX2time,bSymTime);

W = bsxfun(@times,Wspc,reshape(Wtime,[1 1 1 length(Wtime)]));

if bPLOT
    %%
   figure;
   for i = 1:size(W,4), 
      % FULL X,Y,T MOVIE
       subplot(1,2,1);
       imagesc(squeeze(W(:,:,1,i))); 
       Fig.format([],[],[i]); 
       caxis([0 1]);  axis square; 

       % SLICE X,T MOVIE
       subplot(1,2,2);
       ind = floor(size(W,2)/2 + 1);
       plot(squeeze(W(:,:,1,i))); 
       Fig.format(['X-position'],[],['Time=' num2str(i)]); 
       xlim([0 size(W,2)+1]); ylim([0 1.1]); axis square; 
       pause(.1);  
   end
end