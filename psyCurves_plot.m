classdef psyCurves_plot < handle
properties(Access=private)
    bSame
    bLookup
    bSUBJS


end
methods(Access=protected)
    function xtitl=get_xtitl(obj,opts)
        xtitl=obj.EXPS.get_xtitl();
    end
    function ytitl=get_ytitl(obj,opts)
        ytitl='Proportion Cmp. Chosen';
    end
    function titl=get_titl(obj,opts)
        titl=obj.EXPS.get_titl;
    end
    function ititl=get_ititl(obj)
        ititl=cell(obj.m,1);
        for i = 1:obj.m
            ititl{i}=[ 'T = ' num2str(obj.curs{i}.tFit)];
        end
    end
    function [stdX,stdXunq]=get_std(obj,opts)
        stdX=zeros(obj.m,1);
        for i = 1:obj.m
            stdX(i)=unique(obj.curs{i}.stdX);
        end
        stdXunq=unique(stdX);
    end
    function [cmpX,cmpXunq]=get_cmp(obj,opts)
        for i = 1:obj.m
            cmpX(i,:)=unique(obj.curs{i}.cmpX);
        end
        cmpXunq=unique(cmpX);
    end
    function xticks=get_xticks(obj,cmpX)
        % XXX
        %xticks=[-15:2:0];

        %xticks=unique(round(obj.tX,2));

        xticks=unique(cmpX,'rows');
        xticks=cell(size(xticks,1),1);
        for i = 1:size(xticks,1)
            xticks{i}=xticks(i,:);
        end
    end
    function r=get_r(obj,stdXunq,SUBJS)
        if obj.bLookup
            r=numel(stdXunq);
        elseif ~obj.bSUBJS
            r=obj.nSubj;
        elseif iscell(SUBJS) || isnumeric(SUBJS)
            r=length(SUBJS);
        else
            r=[];
        end
    end
    function c=get_c(obj,lookup,stdNums)
        if obj.bSame==1
            c=1;
        elseif obj.bLookup
            s=num2cell(stdNums);
            inds=lookup.lvl('stdInd',s{:}).ret();
            c=numel(unique(inds(:,3)));
            %c=obj.get_lookup_c(lookup);
            %ctitl=obj.get_lookup_ctitle(lookup);
        elseif ~obj.bLookup & obj.bSame==0
            c=obj.m;
        else
            c=[];
        end
    end
    function [rtitl,r]=get_rtitl(obj,SUBJS,stdNums,lookup)
        if obj.bLookup

            if obj.bSame
                rtitl=[];
            else
                s=num2cell(stdNums);
                inds=lookup.lvl('stdInd',s{:}).ret();
                bins=unique(inds(:,3));
                rtitl=arrayfun(@(x) ['Bin ' num2str(x)],bins,'UniformOutput',false);
            end
        elseif ~obj.bSUBJS
               ;
        elseif iscell(SUBJS) || isnumeric(SUBJS)
            rtitl=rtitl(SUBJS);
        else
            rtitl=obj.EXPS.subjs;
        end
    end
    function [ctitl]=get_ctitl(obj,stdXunq)
        if obj.bSame==1
            ctitl=[];
        elseif obj.bLookup
            ctitl=arrayfun(@(x) ['Std ' num2str(x)],stdXunq,'UniformOutput',false);
        elseif ~bLookup & obj.bSame==0
            ctitl=[];
        end
    end
    function lookup=get_lookup(obj,alias)
        lookup=[];
        if obj.bLookup
            lookup=Blk.load_lookup(alias);
        end
    end
    function Opts=get_psy_plot_opts(obj,opts,r,cmpX)
        Opts=struct();
        Opts.xtickFormat=obj.psyXtickFormat;
        Opts.ytickFormat=obj.psyYtickFormat;
        Opts.xtickAngle =obj.psyXtickAngle;
        Opts.ylim       =obj.psyYlim;
        %Opts.xlim       =[];
        %Opts.ax         ='auto';
        Opts.ax='square';
        Opts.position=[3 3 350*r+200 1060 ];
        %Opts.position=[3 3 1060 350*r+200];
        Opts.yticks=[0 .25 .5 .75 1];
        Opts.xtickFormat='%i';
        Opts.xlimSpace=0;
        Opts.xtickAngle=0;
        Opts.bHideXtickLabels=1;
        Opts.bHideXticks=0;
        if obj.bLookup
            Opts.xlimRCBN='C';
        end
        %Opts.xticks=obj.get_xticks(cmpX);
        Opts=Struct.combinePref(opts,Opts);

    end
    function obj = plot_curve_all_p(obj,SUBJS,opts,bSame,alias,stdNums)
         % XXX NOTE
        edges=[.01 .02; .1 .2];

        if ~exist('opts','var') || isempty(opts)
            opts=struct();
        end
        if ~exist('bSame','var') || isempty(bSame)
            bSame=1;
        end
        obj.bSame=bSame;
        if ~exist('alias','var') || isempty(alias)
            bins=[];
            alias=[];
            bLookup=0;
        else
            bLookup=1;
        end
        obj.bLookup=bLookup;
        lookup=obj.get_lookup(alias);

        obj.bSUBJS= exist('SUBJS','var') && ~isempty(SUBJS);

        obj.generate_colors();
        obj.parse_opts();

        xtitl=obj.get_xtitl();
        ytitl=obj.get_ytitl();
        titl=obj.get_titl();
        syms='oooo'; % XXX

        [stdX,stdXunq]=obj.get_std();
        [cmpX,cmpXunq]=obj.get_cmp();

        c=obj.get_c(lookup,stdNums);
        r=obj.get_r(stdXunq,SUBJS);
        Opts=obj.get_psy_plot_opts(opts,r,cmpX);
        Opts.ititl=obj.get_ititl();

        % R AND RTITL
        [rtitl]=obj.get_rtitl(SUBJS,stdNums,lookup);
        [ctitl]=obj.get_ctitl(stdXunq);


        if obj.bSame ~= -1
            sp=SubPlots([c r],xtitl,ytitl,titl,rtitl,ctitl,Opts);
        else
            sp=[];
        end

        %Ymin=[];
        %Ymax=[];
        %Xmin=[];
        %Xmax=[];
        for ind = 1:obj.m
            [~,subj]=ind2sub([obj.m obj.nSubj],ind);
            std=stdX(ind);
            if obj.bSUBJS && ~ismember(subj,SUBJS)
                continue
            elseif obj.bSUBJS
                subj=find(SUBJS==subj);
            end

            % select data
            obj.ind=ind;
            obj.select_cur();
            if isempty(obj.cur)
                continue
            end

            obj.get_col_opts();

            n=obj.get_n(std,stdXunq);
            m=obj.get_m(ind,lookup,stdNums);

            obj.select(sp,m,n);
            obj.plot(syms,subj);

        end
        if obj.bSame ~=-1
            sp.finalize();
        end

    end
    function obj=plot(obj,syms,subj)
        hold on
        obj.plot_curve_p([],syms(subj));
        if obj.bSame == -1 && bLookup
            obj.format_lookup_each(obj,std,Opts);
        end
        hold off
    end
    function obj=select(obj,sp,m,n)
        if obj.bSame ~= -1
            sp.select(m,n);
        else
            Fig.new()
        end
    end
    function m=get_m(obj,ind,lookup,stdNums)
        if obj.bSame ==1
            m=1
        elseif obj.bLookup
            s=num2cell(stdNums);
            inds=lookup.lvl('stdInd',s{:}).ret();
            bin=inds(ind,3);
            bins=unique(inds(:,3));
            m=find(bins==bin);

            %lookup.lvl
            %dk
            %[m,~]=ind2sub([c r],ind);
            %inds=lookup.lvl('stdInd',std).ret();
            %n=inds(3);
        elseif obj.bSame ==0
            m=ind;
        end
    end
    function n=get_n(obj,std,stdXunq)
        if obj.bLookup
            n=find(std==stdXunq);
        else
            n=subj;
        end
    end
    function r=get_lookup_r(obj,lookup)
        tbl=lookup.lvl.ret();
        r=max((tbl(:,2)));
    end
    function c=get_lookup_c(obj,lookup)
        tbl=lookup.lvl.ret();
        c=max((tbl(:,3)));
    end
    function rtitl=get_lookup_rtitle(obj,lookup,stdNums)
        rtitl=cell(obj.m,1);
        for i = 1:obj.m
            inds=lookup.lvl('stdInd',stdNums{i}).ret();
            inds=inds(2);
            rtitl{i}=num2str(obj.EXPS.stdXunqAll{1}(inds));
        end
    end
    function ctitl=get_lookup_ctitle(obj,lookup)
        ctitl=cell(obj.m,1);
        for i = 1:obj.m
            [std,subj]=ind2sub([obj.nStd obj.nSubj],i);
            inds=lookup.lvl('stdInd',std).ret();
            inds=inds(3);
            ctitl{i}=num2str(obj.EXPS.stdXunqAll{2}(inds));
        end
    end

    function format_lookup_each(obj,std,lookup,Opts)
        inds=lookup.lvl('stdInd',std).ret();
        inds=inds(2:end);
        keys=lookup.lvl.KEY;
        keys=keys(2:end);
        titl=[];
        for i = 1:length(obj.EXPS.stdXunqAll)
            titl=[titl ' ' keys{i} ' ' num2str(obj.EXPS.stdXunqAll{i}(inds(i))) newline];
        end
        titl=titl(1:end-1);

        Fig.format(xtitl,ytitl,titl);
        yticks(Opts.yticks);
        xticks(Opts.xticks);
        ylim([0 1]);
        axis square;
    end
    function obj = plot_curve_p(obj,ind,sym)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
            obj.select.cur();
        end
        obj.cur.plot_boot_curve();
        obj.cur.bExtendCurve=0.15;
        obj.cur.plot_data(sym,obj.Opts.color);
        hold off
    end
    function obj=generate_colors(obj)
        if ischar(obj.colormap)
            cmp=eval(obj.colormap);
        elseif isnumeric(obj.colormap)
            cmp=obj.colormap;
        end
        inds=round(linspace(1,size(cmp,1),obj.nStd));
        colors=cmp(inds,:);
        obj.lcolors=colors;
        obj.colors=colors;
    end

end
end
