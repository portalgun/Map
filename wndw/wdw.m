classdef wdw
methods(Static=true)
    function obj=cos(PszRCT,rmpDm,dskDm,symInd)
        if ~exist('symInd','var')
            symInd=[];
        end
        obj=cosWndw(PszRCT,rmpDm,dskDm,symInd)
    end
end
end
