function W = crcWindowXYT(PszXYT,radius,rmpDmT,dskDmT,bPLOT)

% function W = crcWindowXYT(PszXYT,radius,rmpDmT,dskDmT,bPLOT)
%
% example call: W = crcWindowXYT([32 32 48],16,24,24,1);
%
% creates circular window in space with a cosine window (flat top) in time
%
% inputs:
%         PszXYT: size of window in XYT
%         radius: radius of circular window
%         rmpDmT: length of ramp for window in time
%         dskDmT: length of disk for window in time
%         bPLOT : plot or not? 

% MAKE CIRCULAR WINDOW
Wxy = Msk.circle(PszXYT(1:2),radius,[0 0],0);
% CREATE T WINDOW
Wt = cosWindowFlattop([1 PszXYT(3)],dskDmT,rmpDmT,1);
W = bsxfun(@times,Wxy,reshape(Wt,[1 1 length(Wt)]));

if bPLOT==1
    %%
   for i = 1:size(W,3)
       figure;
       clf;
       imagesc(W(:,:,i));
       axis square
       axis xy
       Fig.format(['WX'],['WY'],['RmpT=' num2str(rmpDmT) ', DskT=' num2str(dskDmT)]);
       caxis([0 1])

       suptitle(['t=' num2str(i)],20)
       pause(.05); 
   end
end 

end