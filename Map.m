classdef Map < handle
%TODO reall GAMMA
%TODO DOWNSAMPLE
% XXX interpolate between pixels
properties
    PszXY
    PszRC
    img=cell(1,2)
    LorR
    index

    W=cell(1,2)
    winTex=cell(1,2)
    imgOrig=cell(1,2)
    IbOld

    RMS
    DC
    RMSbino
    DCbino
    RMSmono
    DCmono

    figNum
    fig
    partition %normalizing factor

    dnk=1
    % STATE
    bUpdated
    bAutoUpdate=1
    bContrastFixed=0
    bLuminanceFixed=0
    bWindowed=0
    bFlattened=0
    bNormalized=0
    bCropped=0

    bContrastImg=0
    bEdge=0
    bGamma=0
    monoORbino='bino' %whether to compute stats as single image or independently
    %    Opts.contrast
    %      bByImage %0
    %      monoORbino
    %      DCfix    % .4
    %      RMSfix   %.14
    %      W
    %      nChnl %1

    % OPTSION
    rmsFix
    dcFix
    dnkFix
    bWindow=1
    bFlat=0
    wdwInfo
    flatAnchor=''

    monoORbinoContrast
    monoORbinoFix
end
properties(Hidden=true)
    imgOld=cell(1,2);
    rmsOld
    dcOld
end
methods
    function obj=Map(Limg,Rimg,LorR,WL,WR,index,winTexL,winTexR)
        if nargin < 1
            return
        end

        if nargin >= 1
            if isa(Limg,'ptch')
                obj=init_from_ptch(Limg);
                return
            end
            obj.img{1}=Limg;
            obj.PszXY=size(obj.img{1});
            obj.PszRC=flip(obj.PszXY,2);
        end
        if nargin >= 2
            obj.img{2}=Rimg;
            obj.PszXY=size(obj.img{2});
            obj.PszRC=flip(obj.PszXY,2);
        end
        if nargin >= 3
            obj.LorR=LorR;
        end
        if nargin >=4
            obj.W{1}=WL;
        else
            obj.W{1}=ones(size(obj.img{2}));
        end
        if nargin >=5
            obj.W{2}=WR;
        else
            obj.W{2}=ones(size(obj.img{1}));
        end
        if nargin >=6
            obj.index=index;
        end
        if nargin >=7
            obj.winTex{1}=winTexL;
        else
            obj.winTex{1}=ones(size(obj.img{1}));
        end
        if nargin >=8
            obj.winTex{2}=winTexR;
        else
            obj.winTex{2}=ones(size(obj.img{1}));
        end

        if ~isempty(obj.img{1}) && ~isempty(obj.img{2}) && ~isa(obj,'msk')
            obj=obj.update;
        end
    end
    function obj=init2(obj,bForce,W)
        % NOTE FLATTEN IS HANLDED INDEPENDENTLY
        if nargin < 2 || isempty(bForce)
            bForce=0;
        end
        if nargin >= 3 && ~isempty(W)
            obj.W=W;
        else
            if bForce || (~isempty(obj.wdwInfo) && numel(fieldnames(obj.wdwInfo))~=0)
                obj.gen_window();
            end
        end
        obj.contrast();
        if obj.bWindow
            obj.window();
        end
        if isempty(obj.bFlattened)
            obj.bFlattened=false;
        end

        if obj.bFlat && ~obj.bFlattened
            obj.flatten();
        end
        obj.fix_contrast_bi(obj.rmsFix, obj.dcFix);
        % TODO DNK
        % TODO wdw? obj.
    end
%% FROM PTCH
    % XXX
%% CONVERSION
    function [Limg, Rimg]=split_fun(obj,img)
        Limg=img(:,1:obj.PszRC(1));
        Rimg=img(:,obj.PszRC(1)+1:end);
    end
    function name=name_fun(obj)
        switch class(obj)
            case 'lum'
                name='Luminance';
            case 'edg'
                name='Edge';
            case 'cntrst'
                name='Contrast';
            case 'rnge'
                name='Range';
            otherwise
                name=[];
        end
    end
%% flatten
    function obj=flatten(obj,flatAnchor)
        old=obj.flatAnchor;
        if nargin < 2 && ~isempty(flatAnchor)
            obj.flatAnchor=flatAnchor;
        end
        if obj.bFlattened && isequal(old,obj.flatAnchor)
            return
        elseif obj.bFlattened
            obj.unflatten();
        end

        if isnumeric(obj.flatAnchor) && obj.flatAnchor==0
            return
        elseif isnumeric(obj.flatAnchor) && obj.flatAnchor==1 && ~isempty(obj.LorR)
            flatAnchor=obj.LorR;
        elseif ~isempty(obj.flatAnchor);
            flatAnchor=obj.flatAnchor;
        elseif isempty(obj.flatAnchor)
            flatAnchor='L';
        end
        if strcmp(flatAnchor,'B') && ~isempty(obj.index)  && mod(obj.index,2)==0
            flatAnchor='R';
        elseif strcmp(flatAnchor,'B') && ~isempty(obj.index)
            flatAnchor='L';
        end

        if iscell(flatAnchor) && length(flatAnchor)==1
            flatAnchor=flatAnchor{1};
        end

        if flatAnchor=='L'
            k=1;
            nk=2;
        else
            k=2;
            nk=1;
        end
        obj.imgOld{nk}=obj.img{nk};
        obj.img{nk}=obj.img{k};
        obj.bFlattened=true;;
    end
    function obj=unflatten(obj)
        if ~obj.bFlatttened
            return
        end
        nk=find(~cellfun(@isempty,obj.imgOld));
        obj.img{nk}=obj.imgOld{nk};
        obj.bFlattened=false;
    end
%% NORMALIZE
    function [obj]=normalize(obj)
        if obj.bWindowed
            indL=logical(obj.W{1});
            indR=logical(obj.W{2});
        else
            indL=logical(ones(size(obj.img{1})));
            indR=logical(ones(size(obj.img{2})));
        end
        if isempty(obj.img{2})
            tmp=obj.img{1};
            tmp(~indL)=nan;
            obj.partition=sqrt(nansum(tmp(:).^2));
            obj.img{1}=obj.img{1}/obj.partition;
        elseif isequal(obj.monoORbino,'mono')
            Ltmp=obj.img{1};
            Rtmp=obj.img{2};
            Ltmp(~indL)=nan;
            Rtmp(~indR)=nan;
            obj.partition(1)=sqrt(nansum(Ltmp(:).^2));
            obj.partition(2)=sqrt(nansum(Rtmp(:).^2));
            obj.img{1}=obj.img{1}/obj.partition(1);
            obj.img{2}=obj.img{2}/obj.partition(2);
        elseif isequal(obj.monoORbino,'bino')
            Ltmp=obj.img{1};
            Rtmp=obj.img{2};
            Ltmp(~indL)=nan;
            Rtmp(~indR)=nan;
            tmp=[Ltmp Rtmp];
            obj.partition=sqrt(nansum(tmp(:).^2));
            img=[obj.img{1} obj.img{2}]/obj.partition;
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end
        obj.bNormalized=1;
        if obj.bAutoUpdate
            obj=obj.update;
        end
    end
%% CONTRAST
    function obj=contrast(obj,monoORbino)
        if nargin >= 2 && ~isempty(monoORbino);
           ;
        elseif ~isempty(obj.monoORbinoContrast)
            monoORbino=obj.monoORbinoContrast;
        else
            monoORbino=obj.monoORbino;
        end

        if isempty(obj.img{2})
            [obj.img{1},Ib]=obj.contrast_helper(obj.img{1},obj.W{1});
        elseif strcmp(monoORbino,'L')
            [obj.img{1},Ib(1)]=obj.contrast_helper(obj.img{1},obj.W{1});
            [obj.img{2},Ib(2)]=obj.contrast_helper(obj.img{2},obj.W{2},Ib(1));
        elseif strcmp(monoORbino,'R')
            [obj.img{2},Ib(1)]=obj.contrast_helper(obj.img{2},obj.W{2});
            [obj.img{1},Ib(2)]=obj.contrast_helper(obj.img{1},obj.W{1},Ib(1));
        elseif strcmp(monoORbino,'mono')
            [obj.img{1},Ib(1)]=obj.contrast_helper(obj.img{1},obj.W{1});
            [obj.img{2},Ib(2)]=obj.contrast_helper(obj.img{2},obj.W{2});
        elseif strcmp(monoORbino,'bino')
            [img,Ib]=obj.contrast_helper([obj.img{1} obj.img{2}],[obj.W{1} obj.W{2}]);
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end
        obj.IbOld=Ib;
        obj.bContrastImg=1;
    end
    function [img,Ib]=contrast_helper(obj,img,W,Ib)
        if nargin < 4 || isempty(Ib)
            Ib=sum(img(:).*W(:)./sum(W(:)));
        else
            error('todo')
        end
        img=(img-Ib)/Ib;
    end
    function obj=uncontrast(obj)
        if strcmp(obj.monoORbino,'mono');
            obj.img{1}=obj.uncontrast_helper(obj.img{1});
            obj.img{2}=obj.uncontrast_helper(obj.img{2});
        else strcmp(obj.monoORbino,'bino');
            img=obj.uncontrast_helper([obj.img{1} obj.img{2}]);
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end
        obj.bContrastImg=0;
    end
    function img=uncontrast_helper(obj,img)
        img=img*obj.IbOld+obj.IbOld;
    end
    function obj=unfix_contrast(obj)
        if isequal(obj.monoORbino,'mono')
            error('not implemented')
        end
        obj.fix_contrast_bi(obj.rmsOld{1},[]);
        obj.bContrastFixed=0;
    end
    function obj=unfix_dc(obj)
        if isequal(obj.monoORbino,'mono')
            error('not implemented')
        end
        obj.fix_contrast_bi([],obj.dcOld{1});
        obj.bLuminanceFixed=0;
    end
    function obj=fix_contrast_bi(obj,RMSfix,DCfix,nChnl,monoORbino)
        if nargin >= 5 && ~isempty(monoORbino)
           ;
        elseif ~isempty(obj.monoORbinoFix)
            monoORbino=obj.monoORbinoFix;
        else
            monoORbino=obj.monoORbino;
        end
        if ~obj.bContrastFixed
            bSaveOldRms=1;
        else
            bSaveOldRms=0;

        end
        if ~obj.bLuminanceFixed
            bSaveOldDC=1;
        else
            bSaveOldDC=0;
        end
        if (nargin >=2 && ~isempty(RMSfix))
            obj.bContrastFixed=1;
        else
            RMSfix=[];
        end
        if (nargin >=3  && ~isempty(DCfix)) || obj.bLuminanceFixed
            obj.bLuminanceFixed=1;
        else
            DCfix=[];
        end
        if isempty(RMSfix) && isempty(DCfix)
            return
        end

        if nargin < 4
            nChnl=1;
        end
        dcOld=cell(1,2);
        rmsOld=cell(1,2);
        if isempty(obj.img{2})
            [obj.img{1},rmsOld{1},dcOld{2}]=Map.fix_contrast(obj.img{1},RMSfix,obj.W{1},DCfix,nChnl,obj.bContrastImg,obj.bWindowed);
        elseif isequal(monoORbino,'mono')

            [obj.img{1},rmsOld{1},dcOld{1}]=Map.fix_contrast(obj.img{1},RMSfix,obj.W{1},DCfix,nChnl,obj.bContrastImg,obj.bWindowed);
            [obj.img{2},rmsOld{2},dcOld{2}]=Map.fix_contrast(obj.img{2},RMSfix,obj.W{2},DCfix,nChnl,obj.bContrastImg,obj.bWindowed);
        elseif isequal(monoORbino,'bino')
            [img,rmsOld{1},dcOld{1}]=Map.fix_contrast([obj.img{1} obj.img{2}],RMSfix,[obj.W{1} obj.W{2}],DCfix,nChnl,obj.bContrastImg,obj.bWindowed);
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end
        if DCfix ~= 0
            obj.bContrastImg=0;
        end
        if obj.bAutoUpdate
            obj=obj.update;
        end
        if bSaveOldRms
            obj.rmsOld=rmsOld;
        end
        if bSaveOldDC
            obj.dcOld=dcOld;
        end
    end
%% WINDOWW
    function obj=init_window(obj)
        if isempty(obj.W{1})
            obj.W{1}=ones(size(obj.img{2}));
        end
        if isempty(obj.W{2})
            obj.W{2}=ones(size(obj.img{1}));
        end
    end
    function obj=window(obj)
        if isempty(obj.W{1}) && isempty(obj.W{2})
            error('No window defined');
        end

        if ~obj.bContrastImg
            obj.contrast();
            obj.window_helper;
            obj.uncontrast;
        else
            obj.window_helper;
        end
        obj=obj.update;
    end
    function obj=window_helper(obj)
        % TODO
        W=Map(obj.winTex{1},obj.winTex{2});
        W=W.contrast();
        obj.imgOrig{1}=obj.img{1};
        obj.imgOrig{2}=obj.img{2};
        obj.img{1}=obj.img{1}.*obj.W{1}+(~obj.W{1}.*W.img{1});
        obj.img{2}=obj.img{2}.*obj.W{2}+(~obj.W{2}.*W.img{2});
        obj.bWindowed=1;
    end
    function obj=gen_window(obj)
        flds=fieldnames(obj.wdwInfo);
        w=obj.wdwInfo;
        for i=1:length(flds)
            fld=flds{i};
            val=w.(fld);
            if ~iscell(val)
                w.(fld)={val, val};
                val=w.(fld);
            end
            for k = 1:2
                if ischar(val{k}) && startsWith(val{k},'@')
                    val{k}=val{k}(2:end);
                    %if isprop(obj,val{k})
                        w.(fld){k}=obj.(val{k});
                    %end
                end
            end
        end
        switch obj.wdwInfo.type
        case {'COS','cos'}
            obj.W{1}=cosWdw.new(w.PszRCT{1},w.rmpDm{1},w.dskDm{1},w.symInd{1}); % SLOW
            obj.W{2}=cosWdw.new(w.PszRCT{2},w.rmpDm{2},w.dskDm{2},w.symInd{2});  % SLOW
        otherwise
            error(['Unhandled window type ' obj.winInfo.type]);
        end
    end
    function obj=cos_window(obj,WszRCT,dskDmRCT,rmpDmRCT)
        % TODO
        if numel(WszRCT) == 3
        elseif numel(WszRCT) == 2
            obj.W{1}=cosWindowXY(fliplr(WszRCT),fliplr(dskDmRCT),fliplr(rmpDmRCT),1);
            obj.W{2}=cosWindowXY(fliplr(WszRCT),fliplr(dskDmRCT),fliplr(rmpDmRCT),1);
        end
    end
    function obj=unwindow(obj)
        if ~obj.bWindowed
            return
        end
        obj.img{1}=obj.imgOrig{1};
        obj.img{2}=obj.imgOrig{2};
        obj.imgOrig{1}=[];
        obj.imgOrig{2}=[];
        obj.bWindowed=0;
    end
%% AVERAGE
    function obj=average_vert(obj)
        obj.img{1}=repmat(mean(obj.img{1},1),size(obj.img{1},1),1);
        obj.img{2}=repmat(mean(obj.img{2},1),size(obj.img{2},1),1);
        obj=obj.update;
    end
    function obj=average_hori(obj)
        obj.img{1}=repmat(mean(obj.img{1},1),size(obj.img{1},1),1);
        obj.img{2}=repmat(mean(obj.img{2},1),size(obj.img{2},1),1);
    end
%% CROP
    function obj=crop(obj,PctrRC,PszXY)
        obj.img{1}=cropImgCtrIntrp(obj.img{1},PctrRC,PszXY);
        obj.img{2}=cropImgCtrIntrp(obj.img{2},PctrRC,PszXY);
        obj.bCropped=1;
        obj.PszRC=size(obj.img{1});
        obj.PszXY=fliplr(obj.PszRC);
    end
    function crop_fun(obj,img)
    end
%% EDGE
    function obj = get_edge_strength(obj,tapNum,direction,order)
        if obj.bEdge
            return
        end
        if ~obj.bContrastImg
            obj.contrast;
        end
        if ~obj.bNormalized
            obj.norm;
        end
        if ~exist('tapNum','var')
            tapNum=5;
        end
        if ~exist('direction','var')
            direction='h';
        end
        if ~exist('order','var')
            order=1;
        end
        obj.img{1}=IderivTap(obj.img{1},tapNum,order,direction);
        obj.img{2}=IderivTap(obj.img{2},tapNum,order,direction);
        obj.bEdge=1;
        if obj.bAutoUpdate
            obj.update();
        end
    end
    function obj = get_transition_region(obj,width)
        obj.img{1}=Msk.extremeEdge(obj.img{1},'R',0);
        obj.img{2}=Msk.extremeEdge(obj.img{2},'L',0);

        obj.img{1}=Msk.outline(obj.img{1},width,0);
        obj.img{2}=Msk.outline(obj.img{2},width,0);
    end
%% GAMMA
    function obj=gamma_correct_bi(obj)
        if obj.bGamma
            return
        end
        if isempty(obj.img{2})
            obj.img{1}=Map.gamma_correct(obj.img{1});
        elseif strcmp(obj.monoORbino,'mono')
            obj.img{1}=Map.gamma_correct(obj.img{1});
            obj.img{2}=Map.gamma_correct(obj.img{2});
        elseif strcmp(obj.monoORbino,'bino')
            img=Map.gamma_correct([obj.img{1} obj.img{2}]);
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end
        obj.bGamma=1;
    end
    function obj=gamma_uncorrect_bi(obj)
        if ~obj.bGamma
            return
        end
        if isempty(obj.img{2})
            obj.img{1}=Map.gamma_uncorrect(obj.img{1});
        elseif strcmp(obj.monoORbino,'mono')
            obj.img{1}=Map.gamma_uncorrect(obj.img{1});
            obj.img{2}=Map.gamma_uncorrect(obj.img{2});
        elseif strcmp(obj.monoORbino,'bino')
            img=Map.gamma_uncorrect([obj.img{1} obj.img{2}]);
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end

        obj.bGamma=0;
    end
    function obj=faux_gamma_correct_bi(obj)
        if obj.bGamma
            return
        end
        if isempty(obj.img{2})
            obj.img{1}=Map.faux_gamma_correct(obj.img{1});
        elseif strcmp(obj.monoORbino,'mono')
            obj.img{1}=Map.faux_gamma_correct(obj.img{1});
            obj.img{2}=Map.faux_gamma_correct(obj.img{2});
        elseif strcmp(obj.monoORbino,'bino')
            img=Map.faux_gamma_correct([obj.img{1} obj.img{2}]);
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end
        obj.bGamma=1;
    end
    function obj=faux_gamma_uncorrect_bi(obj)
        if ~obj.bGamma
            return
        end
        if isempty(obj.img{2})
            obj.img{1}=Map.faux_gamma_uncorrect(obj.img{1});
        elseif strcmp(obj.monoORbino,'mono')
            obj.img{1}=Map.faux_gamma_uncorrect(obj.img{1});
            obj.img{2}=Map.faux_gamma_uncorrect(obj.img{2});
        elseif strcmp(obj.monoORbino,'bino')
            img=Map.faux_gamma_uncorrect([obj.img{1} obj.img{2}]);
            [obj.img{1},obj.img{2}]=obj.split_fun(img);
        end

        obj.bGamma=0;
    end
%% PROPS
    function obj=update(obj)
        obj.init_window;
        [obj.RMSmono,obj.DCmono]=obj.rms('mono');
        [obj.RMSbino,obj.DCbino]=obj.rms('bino');
        if strcmp(obj.monoORbino,'bino')
            obj.RMS=obj.RMSbino;
            obj.DC=obj.DCbino;
        elseif strcmp(obj.monoORbino,'mono')
            obj.RMS=obj.RMSmono;
            obj.DC=obj.DCmono;
        end
    end
    % RMS & DC
    function [RMS,DC]=rms(obj,monoORbino)
        if nargin < 2 || isempty(monoORbino)
            monoORbino=obj.monoORbino;
        end

        if isempty(obj.img{2})
            [RMS,DC]=obj.rms_helper(obj.img{1},obj.W{1});
        elseif strcmp(monoORbino,'mono')
            [RMS(1),DC(1)]=obj.rms_helper(obj.img{1},obj.W{1});
            [RMS(2),DC(2)]=obj.rms_helper(obj.img{2},obj.W{2});
        elseif strcmp(monoORbino,'bino')
            [RMS,DC]=obj.rms_helper([obj.img{1},obj.img{2}],[obj.W{1} obj.W{2}]);
        end
    end
    function [RMS,DC]=rms_helper(obj,img,W)
        if (obj.bContrastImg || (obj.bContrastFixed && obj.rmsFix==0)) && obj.bWindowed
            [RMS,DC]=Map.rmsDeviationPreWindowed(img,W);
        elseif obj.bContrastImg || obj.bContrastFixed && obj.rmsFix==0
            [RMS,DC]=Map.rmsDeviation(img,W);
        elseif obj.bWindowed
            [RMS,DC]=Map.rmsContrastPreWindowed(img,W);
        else
            [RMS,DC]=Map.rmsContrast(img,W);
        end
    end
%% PLOT
    function obj=plot(obj)
        if isempty(obj.figNum)
            obj.figNum=Fig.next;
        end
        obj.update();
        obj.fig=figure(obj.figNum);
        imagesc([obj.img{1} obj.img{2}]);
        Fig.format(['RMS ' sprintf('%3.4f',obj.RMS) newline 'DC ' sprintf('%3.4f',obj.DC) ],'',obj.title_fun);
        Fig.formatIm;
        if obj.bEdge
            obj.fig.Colormap=hot;
        end
    end
    function obj= plot_contrast(obj)
        obj=obj.contrast;
        obj.plot;
    end
    function titl=title_fun(obj)
        name=obj.name_fun;
        if ~isempty(name)
            titl=[name ' Image,'];
        else
            titl=[];
        end
        if ~isempty(obj.index)
            titl=[titl ' ' num2str(obj.index)];
        end
        if ~isempty(obj.LorR)
            titl=[titl ' ' obj.LorR ' Anchor'];
        end
    end
    function obj=apply_opts(obj,Opts)
        obj.LorR  =Opts.LorR;
        obj.index =Opts.index;
        obj.rmsFix=Opts.rmsFix;
        obj.dcFix =Opts.dcFix;
        obj.dnkFix=Opts.dnkFix;
        obj.wdwInfo=Opts.wdwInfo;
        obj.bWindow=Opts.bWindow;
        obj.flatAnchor=Opts.flatAnchor;
        obj.monoORbinoContrast=Opts.monoORbinoContrast;
        obj.monoORbinoFix=Opts.monoORbinoFix;
        obj.bFlat=Opts.bFlat;

        obj.imgOld=Opts.imgOld;
        obj.rmsOld=Opts.rmsOld;
        obj.dcOld=Opts.dcOld;
    end
    function Opts=get_opts(obj)
        Opts=struct();
        Opts.bFlat=obj.bFlat;
        Opts.flatAnchor=obj.flatAnchor;
        Opts.LorR=obj.LorR;
        Opts.index=obj.index;
        Opts.rmsFix=obj.rmsFix;
        Opts.dcFix=obj.dcFix;
        Opts.dnkFix=obj.dnkFix;
        Opts.bWindow=obj.bWindow;
        Opts.wdwInfo=obj.wdwInfo;
        Opts.monoORbinoContrast=obj.monoORbinoContrast;
        Opts.monoORbinoFix=obj.monoORbinoFix;

        Opts.imgOld=obj.imgOld;
        Opts.rmsOld=obj.rmsOld;
        Opts.dcOld=obj.dcOld;

    end
end
methods(Static=true)
    function [val1,val2,mn]=getCtrVal(im)
        IszRCD=size(im);
        ctrLow=floor(IszRCD(1:2)./2);
        ctr=ctrLow+1;
        val1=im(ctr(1),ctr(2),:);
        if nargout > 1
            val2=im(ctrLow(1),ctrLow(2),:);
            if nargout > 2
                mn=mean([val1,val2]);
            end
        end
    end
%%- RESIZE
    function map=downsample(map,dnk)
        mapSz=size(map(:,:,1));
        if dnk ~=1
            for i = 1:size(map,3)
                map(:,:,i)=imresize(imresize(map(:,:,i),1/dnk,'bilinear'),mapSz,'bilinear');
            end
        end
    end


%%- CROP
    function [m,exitflag]=crop_f(im,PctrRC,PszRC,interpType)
        if nargin < 4  || isempty(interpType)
            interpType='linear';
        end
        if size(im,3)>1
            m=zeros([PszRC size(im,3)]);
            for i = 1:size(im,3)
                m(:,:,i)=Map.crop_f(im(:,:,i),PctrRC,PszRC,interpType);
            end
            return
        end
        PszXY=flip(PszRC,2);

        switch interpType
        case 'none'
            [m,exitflag] = Map.cropImgCtr(im,round(PctrRC),PszXY,1:size(im,3));
        otherwise
            [m,exitflag]= Map.crop_interp(im,PctrRC,PszXY,interpType);
        end
        if exitflag
            %[size(im) PctrRC PszRC]
            error('Not enough space to crop');
        end
    end

    function [m,exitflag]=crop_interp(Im,PctrRC,PszXY,interpType)
        BszXY   = PszXY + [1+ceil(abs(rem(PctrRC(2),1)))  2 ];
        [m,exitflag]=Map.crop_interp_fun(Im,PctrRC,BszXY,PszXY,interpType);
        if exitflag
            BszXY   = PszXY + [1 1];
            [m,exitflag]=Map.crop_interp_fun(Im,PctrRC,BszXY,PszXY,interpType);
        end
    end
    function [m,exitflag]=crop_interp_fun(Im,PctrRC,BszXY,PszXY,interpType)
        [XbffPix,YbffPix]=meshgrid(1:BszXY(1),1:BszXY(2));

        [crp,exitflag]=Map.cropImgCtr(Im, fix(PctrRC),BszXY);
        if exitflag
            m=0;
            return
        end
        intrp=Map.interp2m( crp, XbffPix+rem(PctrRC(2),1), YbffPix, interpType);
        m=Map.cropImgCtr(intrp,[],PszXY);
    end

    function [crp,exitflag]=cropImgCtr(Im,PctrRC,PszXY,indChnl)
        if nargin < 2 || isempty(PctrRC)
            PctrRC(1) = floor((size(Im,1))/2 + 1);
            PctrRC(2) = floor((size(Im,2))/2 + 1);
        end

        if nargin < 4 || isempty(indChnl)
            indChnl = 1:size(Im,3);
        end
        crd=Map.ctr2crd(PctrRC,PszXY);
        if any(crd <= 0 )
            crp=0;
            exitflag=1;
            return
        else
            exitflag=0;
        end
        crp = Map.cropImg(Im, crd,PszXY,indChnl);
    end
    function crp=cropImg(Im,PcrdRC,PszXY,indChnl)
        crp = Im(PcrdRC(1):(PcrdRC(1)+PszXY(2)-1), ...
                 PcrdRC(2):(PcrdRC(2)+PszXY(1)-1), ...
                 indChnl,:);
    end
    function intrp=interpfun(V,Xq,Yq,interpType)
        % TODO
        for d = 1:size(V,3)
            Vq(:,:,d) = interp2(V(:,:,d),Xq,Yq,interpType);
        end
    end
    function PcrdRC=ctr2crd(PctrRC,PszXY)
        PcrdRC = bsxfun(@minus,PctrRC,floor(flip(PszXY,2)./2));
    end

%%- CONTRAST
    function [IccdRMS,kRMS,DC]=fix_contrast(img,RMSfix,W,DCfix,nChnl,bContrastImg,bWindowed)

        if bWindowed
            DC=mean(img(:));
            %indGd=true(numel(W),1);
            indGd=W(:)>0;
        else
            DC=sum(img(:).*W(:))./sum(W(:));
            indGd=true(numel(W),1);
        end
        %if bContrastImg & bWindowed
        %    Iweb=img;
        if bContrastImg
            %Iweb=img;
            Iweb=img-DC;
        else
            Iweb= (img-DC)./DC;
        end
        if bWindowed
            kRMS = sqrt(sum( (Iweb(indGd).^2)./W(indGd))./sum(W(:)));
        else
            kRMS = sqrt(sum( (Iweb(indGd).^2).*W(indGd))./sum(W(:)));
        end

        if nargin < 2 || isempty(RMSfix)
            RMSfix=kRMS;
        end
        if nargin < 4 || isempty(DCfix)
            DCfix   = DC;
        end

        IccdRMS  = (DCfix.*(RMSfix.*Iweb./kRMS) + DCfix);
    end
    function [RMS,DC] = rmsContrast(I,W)
        % COMPUTE MEAN UNDER WINDOW OF UNWINDOWED IMAGE
        DC   = sum(I(:).*W(:))./sum(W(:));
        % WEBER CONTRAST IMAGE
        Iweb = (I-DC)./DC;
        %Iweb( DC == 0 ) = 0; % when DC==0, meaning all the luminance at the pixels are zero, thus the Iweb, rmsContrast should be zero at the pixels.
        % RMS CONTRAST COMPUTED FROM WEBER CONTRAST IMAGE
        RMS  = sqrt( sum( W(:).*( Iweb(:) ).^2 )./sum(W(:)) );
        DCmin = 0.2;

        if abs(DC) < DCmin
            warning(['rmsContrast: WARNING! mean(input)=' num2str(DC) '. Unstable results if mean(input)<' num2str(DCmin) '!!!. You are probably passing a contrast image. See rmsDeviation.m?'])
        end

    end
    function [RMS,DC]= rmsContrastPreWindowed(Iwdw,W)
        if ~exist('W','var') || isempty(W) W = ones(size(Iwdw)); end
        if ~isequal([size(Iwdw)],size(W))  error(['rmsContrastPreWindowed: size of window W [' num2str(size(W,1)) 'x' num2str(size(W,2)) '] must match size(Iwdw)=[' num2str(size(Iwdw)) ']']); end

        % COMPUTE MEAN
        DC = mean(Iwdw(:));
        %DC = mean(Iwdw(:)./W(:).*sum(W(:)));

        % INDICES W. NON-ZERO WINDOW VALUES
        indGd = W(:)>0;

        % WEBER CONTRAST IMAGE
        IwebWdw = (Iwdw-DC)./DC;

        % RMS CONTRAST COMPUTED FROM WEBER CONTRAST IMAGE
        RMS  = sqrt( sum( (IwebWdw(indGd).^2)./W(indGd) )./sum(W(:)) );
    end

    function [RMS,DC] = rmsDeviation(Izro,W)
        % for zero mean image
        %
        % AVG LUMINANCE UNDER WINDOW
        DC   = sum(Izro(:).*W(:))./sum(W(:));
        % ZERO MEAN IMAGE
        Idev = Izro - DC;
        % RMS DEVIATION COMPUTED FROM MEAN-ZERO IMAGE
        RMS  = sqrt( sum( W(:).*( Idev(:) ).^2 )./sum(W(:)) );
    end
    function [RMS,DC]=rmsDeviationPreWindowed(IzroWdw,W)
        if ~exist('W','var') || isempty(W)   W = ones(size(IzroWdw)); end
        if ~isequal([size(IzroWdw)],size(W)) error(['rmsDeviationPreWindowed: size of window W [' num2str(size(W,1)) 'x' num2str(size(W,2)) '] must match size(IzroWdw)=[' num2str(size(IzroWdw)) ']']); end

        % COMPUTE MEAN
        DC = mean(IzroWdw(:));

        % INDICES W. NON-ZERO WINDOW VALUES
        tol = 1e-4;
        indGd = W(:)>tol;

        %Idev = IzroWdw - DC; % XXX
        Idev=IzroWdw;

        % RMS DEVIATION COMPUTED FROM MEAN-ZERO IMAGE
        %RMS  = sqrt( sum( (( IzroWdw(indGd) ).^2 )./W(indGd))./sum(W(:)) );
        RMS  = sqrt( sum( (( Idev(indGd) ).^2 )./W(indGd))./sum(W(:)) );
    end
    function [Iweb,DC] = contrastImageVec(I,W,bPreWndw)
        % function [Iweb,DC] = contrastImageVec(I,W,bPreWndw)
        %
        %   example calls: Iweb =    contrastImageVec(I)
        %                  Iweb =    contrastImageVec(I,W)
        %                  Iweb = W.*contrastImageVec(I,W)
        %
        % matrix of weber contrast images in column vector form
        % from matrix of intensity images in column vector form
        %
        % I:         matrix of images in column vector form      [ Nd x nStm ]
        % W:         window under which to compute contrast      [ Nd x  1   ]
        % bPreWndw:  boolean indicating whether images have been pre-windowed
        %            1 -> images have     been pre-windowed
        %            0 -> images have NOT been pre-windwoed
        % % % % % % % % % %
        % Iweb:      weber contrast image                        [ Nd x nStm ]
        % DC:        mean of image                               [ 1  x nStm ]

        % INPUT HANDLING
        if ~exist('W','var')        || isempty(W)
           W = ones(size(I,1),1);
        end
        if ~exist('bPreWndw','var') || isempty(bPreWndw)
           bPreWndw = 0;
        end

        tol = 0.2; % MINIMUNM MEAN THAT DOES NOT THROW AN ERROR
        if sum(mean(I) < tol),         disp( ['contrastImageVec: WARNING! mean luminance is less than ' num2str(tol) ' in ' num2str(sum(mean(I) < tol)) ' images... Enter a luminance image']);   end
        if size(W,2) ~= 1,             error(['contrastImageVec: WARNING! window is not in column vector form: size(W)=[' num2str(size(W,1)) ' ' num2str(size(W,2)) ']']);   end
        if size(W,1) ~= size(I,1),     error(['contrastImageVec: WARNING! window size [' num2str(size(W,1)) ' ' num2str(size(W,2)) '] does not match image size [' num2str(size(I,1)) ' ' num2str(size(I,2)) ']']); end
        if sum(W(:)>1) > 1
            error(['contrastImageVec: WARNING! check input W. Values greater than 1.0: max(W(:))=' num2str(max(W(:)))]);
        end


        if bPreWndw == 0 % UNWINDOWED
            % MEAN: EASY TO READ (DC = sum(I(:).*W(:))./sum(W(:))
            DC   = bsxfun(@rdivide,sum(bsxfun(@times,I    ,W) ),sum(W));
            % CONTRAST IMAGE
            Iweb  = bsxfun(@rdivide,    bsxfun(@minus,I    ,DC),DC   );
        elseif bPreWndw == 1 % PREWINDOWED
            % NON-ZERO INDICES
            indGd = W(:)>0;
            % MEAN OF PREWINDOWED IMAGE
            DC   = mean(I,1);
            % CONTRAST IMAGE
            Iweb  = bsxfun(@rdivide,    bsxfun(@minus,I    ,DC) ,DC   );
        end
    end
    function im=gamma_correct(im)
        error('not implemented');
    end
    function im=gamma_uncorrect(im)
        error('not implemented');
    end
    function im=faux_gamma_correct(im)
        im=im.^.4;
    end
    function im=faux_gamma_uncorrect(im)
        im=im.^2.5;
    end
    function im=ptbNormalize(im,bBino)
        if nargin < 2 || isempty(bBino)
            bBino=false;
        end
        if iscell(im)
            mini=zeros(1,2);
            maxi=zeros(1,2);
            for i=1:2
                mini(i)=min(im{i}(:));
                maxi(i)=max(im{i}(:));
            end
            for i = 1:2
                if bBino
                    mini=min(mini);
                    maxi=max(maxi);
                    im{i}=im{i}-mini;
                    im{i}=im{i}/(maxi-mini);
                else
                    im{i}=im{i}-mini(i);
                    im{i}=im{i}/(maxi(i)-mini(i));
                end
            end
        else
            mini=min(im(:));
            maxi=max(im(:));
            im=im-mini;
            im=im/(maxi-mini);
        end
    end
    function map=mapToMode(moude,map)
        if strcmp(moude,'sng')
            return
        elseif strcmp(moude,'sbs')
            map=im_to_sbs(map);
        elseif strcmp(moude,'ana')
            map{1}=stereoAnaglyph(map{1},map{2});
            map{2}=map{1};
        else
            error(['Unhandled image mode ' moude]);
        end
        function sbs=im_to_sbs(map)
            sbs=cell(1,2);
            sbs{1}=[map{1} map{2}];
            sbs{2}=sbs{1};
        end
    end


%%- INTERP
    function xyzM=interp(PPxy,IppXm,IppYm)
        xyzM=Map.interp2(IppXm,IppYm,X,Y,PPxy(:,1),PPxy(:,2));
    end
    %function PPxy=revInterp(xyzM,IppXm,IppYm,X,Y)
    %    PPxy=Map.revImgInterp2(IppXm,IppYm,X,Y,xyzM(:,1),xyzM(:,2));
    %end
    function out=gridInterpXY(F,VRC)
        out=zeros(size(VRC));
        out(:,1)=F{1}(VRC(:,2));
        out(:,2)=F{2}(VRC(:,1));
    end
    function [Xitp,Yitp]=revInterp2(Xm,Ym,X,Y,Vxm,Vym)
        % XXX CHECK MINUS -0.5
        %ind=size(X,2)/2
        %X(1,ind)
        %X(1,ind+1)
        if Arr.ndim(Vxm)==1
            Xitp=interp1(Xm(1,:),X(1,:)-0.5,Vxm,'linear');
            Yitp=interp1(Ym(:,1),Y(:,1)-0.5,Vym,'linear');
        else
            Xitp=interp2(Xm,Ym,X(:)-0.5,Vxm,Vym,'linear');
            Yitp=interp2(Xm,Ym,Y(:)-0.5,Vxm,Vym,'linear');
        end
    end
    function [Xitp,Yitp]=interp2(X,Y,Xm,Ym,Vx,Vy)

    % X Y in pixels
    % Ym Xm is meters projection plane
    %
        if Arr.ndim(Vx)==1
            Xitp=interp1(X(1,:),Xm(1,:),Vx,'linear');
            Yitp=interp1(Y(:,1),Ym(:,1),Vy,'linear');
        else
            Xitp=interp2(X,Y,Xm(:),Vx,Vy,'linear');
            Yitp=interp2(X,Y,Ym(:),Vx,Vy,'linear');
        end
    end
    function vq = interp1Clip(x,v,xq)
    % Interpolation for edges

        vq = interp1(x,v,xq);

        [XMax, idxVMax] = max(x);
        [XMin, idxVMin] = min(x);

        idxMax = xq > XMax;
        idxMin = xq < XMin;

        vq(idxMax) = v(idxVMax);
        vq(idxMin) = v(idxVMin);
    end
    function Vq = interp2m(varargin)

        % function Vq = interp2m(varargin)
        %
        % wrapper function for interp2 that handles matrix inputs
        % especially useful for interpolating RGB images

        if     nargin == 3
            V  = varargin{1};
            Xq = varargin{2};
            Yq = varargin{3};
            % INTERPOLATED VALUES
            for d = 1:size(V,3)
            Vq(:,:,d) = interp2(V(:,:,d),Xq,Yq);
            end
        elseif nargin == 4
            V  = varargin{1};
            Xq = varargin{2};
            Yq = varargin{3};
            interpType = varargin{4};
            % INTERPOLATED VALUES
            for d = 1:size(V,3)
            Vq(:,:,d) = interp2(V(:,:,d),Xq,Yq,interpType);
            end
        elseif nargin == 5
            X  = varargin{1};
            Y  = varargin{2};
            V  = varargin{3};
            Xq = varargin{4};
            Yq = varargin{5};
            % INTERPOLATED VALUES
            for d = 1:size(V,3)
            Vq(:,:,d) = interp2(X,Y,V(:,:,d),Xq,Yq);
            end
        elseif nargin == 6
            X  = varargin{1};
            Y  = varargin{2};
            V  = varargin{3};
            Xq = varargin{4};
            Yq = varargin{5};
            interpType = varargin{6};
            % INTERPOLATED VALUES
            for d = 1:size(V,3)
            Vq(:,:,d) = interp2(X,Y,V(:,:,d),Xq,Yq,interpType);
            end
        else
        error(['interp2m: WARNING! unhandled number of input arguments: ' num2str(nargin)]);
        end
    end
    function dots=view_dot_ptchs(P,kernSz)
        j
        if nargin < 1
            P=ptchs.getRaw('all');
        end
        fname=[P.get_dir() '_CP_verify_dot_' strrep(Num.toStr(kernSz),',','-') '_.mat'];
        S=load(fname);
        dots=S.dots;

        bInd=true(size(P.fnames));
        %bInd=P.idx.B==1;
        nInd=main.newInd(P);

        bBd=~nInd & bInd &(P.Flags.seen | P.Flags.other) &  P.Flags.bad;
        bGd=~nInd & bInd &(P.Flags.seen | P.Flags.other) & ~P.Flags.bad;

        sum(bGd)
        sum(bBd)
        gdDots=dots(bGd);
        bdDots=dots(bBd);
        sum(gdDots < .99)
        sum(bdDots < .99)
        bins=linspace(.99,1,50);
        %bins=linspace(0,1,50);

        figure(1)

        hist(gdDots,bins);
        title('bGood');
        xlim([bins(1) bins(end)]);

        figure(2)
        hist(bdDots,bins);
        xlim([bins(1) bins(end)]);
        title('bBad');

    end


end
end
