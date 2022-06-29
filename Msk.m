classdef Msk < handle
methods(Static)

    function I = circle(PszXY,radius,PctrXY,bPLOT)

    % function I = Msk.circle(PszXY,radius,PctrXY,bPLOT)
    %
    %   example call: Msk.circle([500 500],200,[0 0],1);
    %
    % create circle with arbitrary radius and center with values 0 and 1
    %
    % PszXY:  patch size in pixels          [1 x 2]
    % radius: radius in pixels
    % PctrXY: circle center in pixels       [1 x 2]
    % bPLOT:  plot or not
    % %%%%%%%%%%%%%%%%%%
    % I:      circle image

        if numel(PszXY)==1
            PszXY=[PszXY PszXY];
        end
        if nargin < 2
            radius=floor(min(PszXY(:))/2);
        end
        if nargin < 2
            PctrXY=floor(PszXY./2);
        end

        if size(PszXY,1) == 1 && size(PszXY,2) == 1
            [X,Y] = meshgrid(Wave.smpPos(1,PszXY));
        elseif size(PszXY,1) == 1 && size(PszXY,2) > 1
            [X,Y] = meshgrid(Wave.smpPos(1,PszXY(1)),Wave.smpPos(1,PszXY(2)));
        end
        R = sqrt((X-PctrXY(1)).^2 + (Y-PctrXY(2)).^2);

        I = zeros(fliplr(PszXY));
        I(R <= radius) = 1;

        return
        if bPLOT
            figure;
            imagesc(I);
            axis image;
            Fig.format([],[],['Radius=' num2str(radius) 'pix']);
            axis xy;
        end
    end

    function prbCrpInd = circle2(PszXY,ctrXY,r,bPlot)
    %     PszXY=[50 50]
    %     ctrXY=[ 10 10; 15 15 ];
    %     r=3;
    %     bPlot=1;
    %     circle2(PszXY,ctrXY,r,bPlot)
    %Returns logical map of circle of radius r in image with dimensions PszXY
        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=0;
        end
        numel(ctrXY);

        if numel(ctrXY) >= 2
            ind=sub2ind(fliplr(PszXY),ctrXY(:,2),ctrXY(:,1));
            im=zeros(fliplr(PszXY));
            im(ind)=1;
            prbCrpInd = imdilate(im, strel('sphere',r));
            %prbCrpInd = imdilate(im, strel('disk',r,r));
            %prbCrpInd = imdilate(im, strel('rectangle',[1 r*2])); %hori
            %prbCrpInd = imdilate(im, strel('rectangle',[r*2 1])); %hori
        else
            [X Y] = meshgrid(1:PszXY(1), 1:PszXY(2));
            prbCrpInd = (Y - ctrXY(:,2)).^2 + (X - ctrXY(:,1)).^2 <= r.^2;
        end
        if bPlot
            imagesc(prbCrpInd);
            xticks(1:PszXY(1));
            yticks(1:PszXY(2));
            grid on;
            %imagesc(out);
        end
    end
    function prbCrpInd = circle2Vec(PszRC,ctrRC,r,bPlot)
    % function prbCrpInd = circle2(PszXY,ctrXY,r,bPlot)
    %
    %Returns logical map of circle of radius r in image with dimensions PszXY
    %
    % example call
    %     PszXY=[50 50];
    %     ctrXY=[ 10 10; 15 15];
    %     r=3;
    %     bPlot=1
    %     circle2(PszXY,ctrXY,r,bPlot)
    %
        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=0;
        end
        if size(ctrRC,2) == 2
            inds=sub2ind(PszRC,ctrRC(:,1),ctrRC(:,2));
        elseif size(ctrRC,2) == 2
            inds=ctrRC;
        end

        prbCrbInd=zeros(PszRC);
        %[X Y] = meshgrid(1:PszXY(1), 1:PszXY(2));
        X(inds)=1;
        Y(inds)=1;

        prbCrpInd = (Y - ctrXY(2)).^2 + (X - ctrXY(1)).^2 <= r.^2;
        if bPlot
            imagesc(prbCrpInd);
        end
    end
    function prbCrpInd = rec(PszXY,ctrXY,h,w,bPlot)
        % Msk.rec([100 100],[20 20; 50 50],4,2,1)
        %Returns logical map of circle of radius r in image with dimensions PszXY
        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=0;
        end

        %prbCrpInd = (Y - ctrXY(2)).^2 + (X - ctrXY(1)).^2 <= r.^2;
        if size(ctrXY,1) >= 2
            ind=sub2ind(fliplr(PszXY),ctrXY(:,2),ctrXY(:,1));
            im=zeros(fliplr(PszXY));
            im(ind)=1;
            prbCrpInd = imdilate(im, strel('rectangle',[h w])); %horie
        else
            [X Y] = meshgrid(1:PszXY(1), 1:PszXY(2));
            prbCrpInd= abs(Y-ctrXY(2)) <= h & abs(X-ctrXY(1)) <= w;
        end

        if bPlot
            imagesc(prbCrpInd);
        end
    end
    function binMsk=cuboid(PsXYZ,ctrXYZ,h,w,d,bPlot)
        %PszXYZ=[10 10 20];
        %A=intersectLineWithGrid([0 0 0],[10 10 20],PszXYZ);
        %B=intersectLineWithGrid([0 0 20],[10 10 0],PszXYZ);
        %ctrXYZ=[A; B];
        %h=4;
        %w=2;
        %d=1;
        %bPlot=1;
        %bInd=logicalCuboid(PszXYZ,ctrXYZ,h,w,d,bPlot)

        ind=sub2ind(PszXYZ,ctrXYZ(:,1),ctrXYZ(:,2),ctrXYZ(:,3));
        PszRCT=transpose(permute(PszXYZ,[2 1 3]));
        im=zeros(PszRCT);
        im(ind)=1;
        bInd=imdilate(im,strel('cuboid',[h w d]));
        if exist('bPlot','var') && isequal(bPlot,1)
            for i = 1:size(bInd,3)
                imagesc(bInd(:,:,i));
                pause(.2)
                drawnow
            end
        end
    end
    function new = outline(BW,thickness,bNoEdge,bPlot)
        if ~exist('bNoEdge','var') || isempty(bNoEdge)
            bNoEdge=0;
        end
        if ~exist('thickness','var') || isempty(thickness)
            thickness=1;
        end
        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=0;
        end

        if bNoEdge
            Bn=[BW(1,:); BW; BW(end,:)];
            Bn=[Bn(:,1), Bn, Bn(:,end)];
        else
            Bn=BW;
        end

        PszRC=size(Bn);
        B=bwboundaries(Bn);
        B=vertcat(B{:});
        new=zeros(PszRC);
        for i = 1:size(B,1)
            new(B(i,1),B(i,2))=1;
        end

        if bNoEdge
            new=new(:,2:end-1);
            new=new(2:end-1,:);
        end

        if thickness>0
            new=Msk.surround(new,thickness,2);
        end
        if bPlot
            plot_fun(new,BW);
        end
    end
    function imgNew = surround(img,expandLength,dim)
        % logical map expanding in dimension dim at length of expandLength from true elements
        %function imgNew = Msk.surround(img,expandLength,dim)
        %
        % example call:
        %           A=zeros(10);
        %           A(25)=1;
        %           A(38)=1;
        %           Msk.surround(A,2,1)
            assert(mod(expandLength,1)==0,'expandLength must be an integer value');
        imgNew=zeros(size(img));
        PszRC=size(img);
        if ~exist('dim','var') || isempty(dim)
            dim=1;
        end

        range=[-expandLength:expandLength]';

        [Ln,Lm]=find(img);
        l=length(Ln);
        Ln=repelem(Ln,length(range),1);
        Lm=repelem(Lm,length(range),1);

        r=repmat(range,l,1);

        if dim==1
            Ln=Ln+r;
            rmInd=(Ln <= 0 | Ln > PszRC(1));
        elseif dim==2
            Lm=Lm+r;
            rmInd=(Lm <= 0 | Lm > PszRC(2));
        end

        Lm(rmInd)=[];
        Ln(rmInd)=[];

        ind=sub2ind(PszRC,Ln,Lm);
        imgNew(ind)=1;
    end

    function new = edge(map,LorR,bPlot)
        if LorR=='L'
            new=diff(map,1,2)>=1;
        elseif LorR=='R'
            new=diff(map,1,2)<=-1;
        end
        new=[zeros(size(map,1),1) new];

        if bPlot
            plot_fun(new,map);
        end
    end
    function new=extremeEdge(map,LorR,bPlot)
    %LorR refers to the edge
        if LorR=='R'
            map=fliplr(map);
        end
        %new=~cumprod(~bwperim(map),2)
        new=diff(~cumprod(~map,2),1,2);
        new=[zeros(size(map,1),1) new];
        if LorR=='R'
            map=fliplr(map);
            new=fliplr(new);
        end

        if bPlot
            plot_fun(new,map);
        end
    end
    function new = shift(map,r,c,bPlot)

        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=0;
        end

        new=map;

        if 0>r
            r=abs(r);
            szR=size(new(1:r,:));
            new(1:r,:)=[];
            new=[new; zeros(szR)];
        elseif 0<r
            szR=size(new(end-r:end,:));
            new(end-r:end,:)=[];
            new=[zeros(szR);new];
        end

        if 0>c
            c=abs(c);
            szR=size(new(:,1:c));
            new(:,1:c)=[];
            new=[new, zeros(szR)];
        elseif 0<c
            szR=size(new(:,end-c:end));
            new(:,end-c:end)=[];
            new=[zeros(szR),new];
        end

        if bPlot
            plot_fun(new,map);
        end
    end
    function out=maxRunLengthFast(A,dim,val,thresh)
        if nargin < 2
            dim=[];
        end
        if nargin < 3
            val=uint32(1);
        end
        if nargin < 4
            thresh=0;
        end

        [vals,cnts,inds]=Msk.run_length_cpp(A,dim);


        ind=(vals==val & cnts > thresh);
        sinds=single(inds);
        subRC=Index.toSub(sinds(ind),size(A));

        out=accumarray(subRC(:,1),cnts(ind),[],@max);
    end
    function out=nRunLengthFast(A,dim,val,thresh)
        if nargin < 2 || isempty(dim)
            dim=[];
        end
        if nargin < 3 || isempty(val)
            val=uint32(1);
        end
        if nargin < 4 || isempty(thresh)
            thresh=0;
        end

        [vals,cnts,inds]=Msk.run_length_cpp(A,dim);

        ind=(vals==val & cnts > thresh);
        sinds=single(inds);
        subRC=Index.toSub(sinds(ind),size(A));
        %subRC=[rem(dinds(ind),size(A,1))+1 inds(ind)./size(A,1)]

        out=accumarray(subRC(:,1),subRC(:,2),[],@numel);
        if isempty(out)
            out=0;
        end
    end
    function [oneCnt]=maxRunLength(A,dim,thresh)
        if nargin < 3
            thresh=1;
            if nargin < 2
                dim=1;
            end
        end
        if sum(A,'all') == prod(size(A))
            oenCnt=size(A,1);
        elseif sum(A,'all') == 0
            oneCnt=0;
            return
        end
        if dim==1
            A=transpose(A);
        end
        runCnt=Msk.getRunCnt(A);
        d=diff([ones(size(A,1),1),A],[],2);


        %ONE
        lindOneChg=((A==1 & d== 1));
        lindOneChg((A(:,1)==1),1)=1;
        oneChg=Msk.getChg(lindOneChg);

        oneCnt=max(cellfun(@(x) sum(x >= thresh),oneChg));
    end
    function [oneCnt,zerCnt,indOneChg,indZerChg]=runLength(A,dim)
        %oneCnt     - length of each count of 1s (cell nx1)
        %zerCnt     - length of each count of 0s (cell nx1)
        %indOneChg - index of one change points (cell nx1)
        %indOneChg - index of zero change points (cell nx1)

        if nargin < 3
            thresh=1;
            if nargin < 2
                dim=1;
            end
        end
        if dim==1
            A=transpose(A);
        end

        PszRC=size(A);
        if sum(A,'all') == prod(PszRC)
            oneCnt=repmat({PszRC(2)},PszRC(1),1);
            zerCnt=repmat({0},PszRC(1),1);
            indOneChg=repmat({1},PszRC(1),1);
            indZerChg=repmat({},PszRC(1),1);
            return
        elseif sum(A,'all') == 0
            zerCnt=repmat({PszRC(2)},PszRC(1),1);
            oneCnt=repmat({0},PszRC(1),1);
            indZerChg=repmat({1},PszRC(1),1);
            indOneChg=repmat({},PszRC(1),1);
            return
        end

        runCnt=Msk.getRunCnt(A);

        d=diff([ones(size(A,1),1),A],[],2);

        %ONE
        lindOneChg=((A==1 & d== 1));
        lindOneChg((A(:,1)==1),1)=1;
        oneChg=Msk.getChg(lindOneChg);

        %ZERO
        lindZerChg=((A==0 & d==-1));

        %BOTH
        lindChg=(lindOneChg | lindZerChg);
        chg   =Msk.getChg(lindChg);

        %oneCnt=cellfun(@(x,y,z) z(ismember(x,y)), chg,indOneChg,runCnt,'UniformOutput',false);
        oneCnt=cellfun(@numel,oneChg,'UniformOutput',false);
        if nargout > 1
            indZerChg=cellfun(@(x,y) x(~ismember(x,y)), chg,oneChg,'UniformOutput',false);
            zerCnt=cellfun(@(x,y,z) z(ismember(x,y)), chg,indZerChg,runCnt,'UniformOutput',false);
            if nargout > 2
                indOneChg=cellfun(@(x,y) x(ismember(x,y)),  chg,oneChg,'UniformOutput',false);
                if nargout > 3
                end
            end
        end
    end
    function [indChg,runVal,runCnt] = runLengthSimple(x)

        % function [indChg,runVal,runCnt] = runLengthEncoder(x)
        %
        %   example call: [indChg,runVal,runCnt] = runLengthEncoder([0 0 1 1 1 0 1 1 0 0])
        %
        % uses run length encoding algorithm to compress string of booleans
        %
        % x:         vector of booleans         [ 1 x n ] or [ n x 1 ]
        %%%%%%%%%%%%%%%%%%%%
        % indChg:    indices of change point    [ 1 x nRun+1 ]
        %            [1 length(x)]
        % runVal:    value   of each run        [ 1 x nRun   ]
        % runCnt:    length  of each run        [ 1 x nRun   ]

        if size(x,1) > size(x,2), x = x'; end % if x is a column vector, tronspose
        if size(x,1) > 1 && size(x,2) > 1, error(['runLengthEncoder: WARNING! x must be a vector. Currently size(x)=[' num2str(size(x)) ']' ]); end

        i = [ find(x(1:end-1) ~= x(2:end)) length(x) ];
        runVal    = x(i);
        runCnt = diff([ 0 i ]);
        indChg    = [1 cumsum(runCnt(1:end-1))+1 length(x)];
    end

    function [new,MaskBG,IMG]=cutAndFill(map,IctrRC,mask,insert,bNormalize,bPlot)
        if ~exist('insert','var') || isempty(insert)
            insert=zeros(size(mask));
        end
        if ~exist('bPlot') || isempty(bPlot)
            bPlot=0;
        end
        if ~isequal(size(insert),size(mask))
            insert=imresize(insert,size(mask),'method','bilinear');
        end
        if ~exist('bNormalize') || isempty(bNormalize)
            bNormalize=1;
        end

        [R,C]=find(mask);
        RR=R+(IctrRC(1)-mean(R))-1;
        CC=C+(IctrRC(2)-mean(C))-1;
        RC=[RR CC];

        if isempty(RC)
            shift=[0 0];
        else
            shift=mod(RC(1,:),1);
        end

        x=0:size(map,2);
        y=0:size(map,1);
        [x,y]=meshgrid(x,y);

        %X=x
        %Y=y
        X=x+shift(2);
        Y=y+shift(1);

        val=~ismember([Y(:),X(:)],RC,'rows');
        if sum(sum(~val))==0
            X=x+shift(2);
            Y=y+shift(1);
            val=~ismember([Y(:),X(:)],RC,'rows');
        end
        val=reshape(val,size(X));

        [IR,IC]=find(~val);
        MaskBG=zeros(size(map));
        MaskBG(IR,IC)=mask(R,C);
        MaskBG=imtranslate(MaskBG,-1*shift,'FillValues',0);
        %MaskBG=MaskBG>.5;
        MaskBGtmp=ones(size(MaskBG));

        [IR,IC]=find(~val);
        IMG=zeros(size(map));
        IMG(IR,IC)=insert(R,C);
        IMG=imtranslate(IMG,-1*shift,'FillValues',0);
        IMG=IMG.*MaskBGtmp;
        %new=map.*(1-MaskBG)+IMG.*MaskBG;

        new=map.*(1-MaskBG)+IMG.*MaskBG;

        if bNormalize
            MaskBG=MaskBG/(max(max(MaskBG)));
            IMG=IMG/max(max(IMG));
        end


        %% Fix nans
        ind=isnan(new);
        new(ind)=map(ind);
        MaskBG(ind)=1;

        if bPlot==1
            plot_fun2(IctrRC,new,map,MaskBG,IMG);
        end

    end
    function plot_fun2(IctrRC,new,map,Mask,IMG)
        figure(22)
        subplot(2,2,1)
        imagesc(map)
        Fig.formatIm;
        colorbar
        title('map')
        c=caxis;

        subplot(2,2,2)
        imagesc(IMG)
        Fig.formatIm;
        colorbar
        title('IMG')
        hold on
        plot(IctrRC(2),IctrRC(1),'.r')
        caxis(c);


        subplot(2,2,3)
        imagesc(new)
        Fig.formatIm;
        colorbar
        title('new')
        hold on
        plot(IctrRC(2),IctrRC(1),'.r')
        caxis(c);

        subplot(2,2,4)
        imagesc(Mask)
        Fig.formatIm;
        colorbar
        title('Mask')
        hold on
        plot(IctrRC(2),IctrRC(1),'.r')
        caxis(c);
    end

    function plot_fun(new,map)
        figure(22)
        subplot(1,2,1)
        imagesc(map);
        Fig.formatIm;
        subplot(1,2,2)
        imagesc(new);
        Fig.formatIm;
    end
end
methods(Static, Access=private)
    function [start,cnt,inds]=run_length_cpp(A,dim)
        if nargin < 2 || isempty(dim)
            dim=2;
        end
        if dim==1
            A=A';
        end
        [start,cnt,inds]=rleMex(uint16(A));

    end
    function chg = getChg(lindChgMap)
        PszRC=size(lindChgMap);
        [M,N]=find(lindChgMap');

        % SLOW 2
        [~,i]=unique(N,'last');
        ind=[[1; i(1:end-1)+1] i];

        % SLOW 1
        chgRaw=arrayfun(@(r,c) M(r:c),ind(:,1),ind(:,2),'UniformOutput',false);

        %FILL IN MISSING
        if length(chgRaw)~=PszRC(1)
            chg=cell(PszRC(1),1);
            lind=transpose(ismember(1:PszRC(1),N));
            chg(lind)=chgRaw;
        else
            chg=chgRaw;
        end
    end
    function  runCnt =getRunCnt(A)
        [M,N]=find([ ones(size(A,1),1) diff(A,[],2) ones(size(A,1),1) ]');
        % SLOW 2
        [~,i]=unique(N,'last');

        ind=[[1; i(1:end-1)+1] i];

        % SLOW 1
        runCnt=arrayfun(@(r,c) diff(M(r:c)),ind(:,1),ind(:,2),'UniformOutput',false);
    end
end
end
