classdef Tx < handle
properties
    self
    genName
    opts
    PszRC
    im
end
methods
    function obj=Tx(PszRC,genName,im);
        if ~exist('im','var')
            im=[];
        end
        obj.im=im;
        obj.PszRC=PszRC;
        obj.genName=genName;
        [obj.self,obj.opts]=Tx.gen(PszRC,genName,im);
    end
    function plot(obj)
        imagesc(obj.self);
        Fig.formatIm();
        caxis([0,1]);
    end
end
methods(Static=true)
    function out=isgen(list)
        %0.5
        %f2.0.5
        out=cellfun(@fun,list);
        function out=fun(in)
            out=isnum(in(2:end));
        end
    end
    function [tex,opts]=gen(PszRC,genName,im)
        % f (exponent)
        % sd
        % rms
        % dc

        opts=parse_str(genName);

        if isnan(opts.f) && isnan(opts.rms)
            tex=ones(PszRC).*opts.dc;
            return
        elseif isnan(opts.f) && ~isnan(opts.rms)
            error('tex: if no exponent specific, cannot set rms');
        elseif ~isnan(opts.f)
            tex=Noise.img(fliplr(PszRC),opts.f,0,opts.sd);
        end

        function opts=parse_str(genName)
            opts=struct();
            opts.f=nan;
            opts.dc=0.5;
            opts.rms=nan;
            opts.sd=[];
            flds=fieldnames(opts);

            spl=strsplit(genName,'_');

            for i = 1:length(spl)
                spl=spl{i};
                opt=Str.RE.match(spl,'[a-z]*');
                val=strrep(spl,opt,'');
                if Str.Num.isReal(val);
                    val=str2double(val);
                end
                if ismember(opt,flds)
                    opts.(opt)=val;
                else
                    error(['unrecognized option ' opt]);
                end

            end
        end

    end
end
end
