classdef PsyCurvePlot < handle
properties
    %colors=[[18,133,158];[49,97,107];[45,209,148];[215,103,100];[158,18,73];[18,133,158]];
    %lcolors=colors ./ 255;
    PsyCurve

    nX % curve resolution
    color % line color
    lineColor
    lineWidth
    CIAlpha

    markersize
    markerface

    XName
    XUnits

    bPlotCI
    bExtendCurve
    bTitle

    Ymin
    Ymax
    Xmin
    Xmax
end
methods(Static)
    function P=getP()
        P={...
            'bPlotCI',[],'isbinary';
            'bTitle',1,'isbinary';
            'lineWidth',2,'Num.is';
            'markersize',12,'Num.is';
            'markerface','w','@(x) true';
            'color','k','@(x) true';
            'lineColor','k','@(x) true';
            'CIAlpha',.4,'Num.is';
            'shape','o','@(x) true';
            'XName','X','ischar';
            'Xunits','','ischar_e';
            'nX',200,'Num.is';
            'bExtendCurve',0,'isbinary';
            'shape',0,'isbinary';
        };
    end
end
methods
    function obj=PsyCurvePlot(parent,varargin)
        if nargin < 1
            return
        end
        obj.PsyCurve=parent;
        if nargin < 2
            return
        end
        obj.parse(varargin{:});
    end
    function parse(obj,varargin)
        Args.parse(obj,obj.getP,varargin{:});
    end
    function plot(obj,meas,units,mult,xfrmt,opts);
        obj.plot_();

        obj.Format(meas,units,mult,xfrmt,opts);
        hold off;
    end
    function format(obj,meas,units,mult,xfrmt,opts)
        obj.xylim([],opts);
        obj.ylabel();
        obj.rylabel(mult,opts);
        obj.cmplabel(meas,units);
        obj.cmpticks(mult,xfrmt);

        Axis.format;
    end
    function [] = plot_all(obj)
        Fig.new();
        hold off;
        subPlot([2,2],1,1);
        obj.plot_boot_curve();
        obj.plot_boot_params(2);
        subPlot([2 2],2,1);
        obj.plot_boot_DP();
    end
    function [] = plot_(obj,opts)
        if nargin < 2
            opts=[];
        end
        obj.plot_boot_curve(opts);
        hold on;
        obj.plot_data([],[],opts);
    end
    function [] = plot_params(obj)
        Fig.new();
        hold off;
        obj.plot_boot_params();
    end
    function [] = plot_DP(obj)
        Fig.new();
        hold off;
        obj.plot_boot_DP();
    end
end
methods(Hidden = true)
    function [] = plot_boot_params(obj,c)
        if ~exist('c','var') || isempty(c)
            c=1;
        end

        obj.get_bP();
        sSz=[sum(obj.bP)+1,c];

        subPlot(sSz,1,c);
        obj.plot_boot_T_p();

        i=1;
        if obj.bP(i)
            subPlot(sSz,i+1,c);
            obj.plot_boot_mu_p();
        end

        i=2;
        if obj.bP(i)
            subPlot(sSz,i+1,c);
            obj.plot_boot_sigma_p();
        end

        i=3;
        if obj.bP(i)
            subPlot(sSz,i+1,c);
            obj.plot_boot_bet_p();
        end


    end

    function [Xdat,Ydat]=get_dataXY(obj)
        Xdat = transpose(unique(obj.cmpX));
        % PROPORTION CMP CHOSEN
        Ydat=zeros(length(Xdat),1);
        for i = 1:length(Xdat)
            ind = obj.cmpX(:) == Xdat(i);
            Ydat(i) = mean(obj.RCmpChs(ind));
        end
    end
    function plot_data(obj,sym,color,opts)
        if nargin >=2 && isstruct(sym) || nargin >=4 && isstruct(opts)
            if nargin >=2 && isstruct(sym)
                opts=sym;
            end
            if isfield(opts,'Marker')
                sym=opts.Marker;
            else
                sym=obj.Shape;
            end
            if isfield(opts,'MarkerColor')
                color=opts.MarkerColor;
            else
                color=obj.lcolor;
            end
            if isfield(opts,'MarkerSize')
                markersize=opts.MarkerSize;
            else
                markersize=obj.markersize;
            end
            if isfield(opts,'MarkerFace')
                markerface=opts.MarkerFace;
            else
                markerface=obj.markerface;
            end
            if isfield(opts,'LineWidth')
                LineWidth=opts.LineWidth;
            else
                LineWidth=obj.LineWidth;
            end
            if isfield(opts,'bAbs')
                bAbs=opts.bAbs;
            else
                bAbs=false;
            end
        else
            if ~exist('sym','var') || isempty(sym)
                sym=obj.shape;
            end
            if ~exist('color','var') || isempty(color)
                color=obj.lcolor;
            end
            markersize=obj.markersize;
            markerface=obj.markerface;
            LineWidth=obj.LineWidth;
            bAbs=false;
        end
        [Xdat,Ydat]=obj.get_dataXY();
        if bAbs && all(Xdat < 0)
            Xdat=abs(Xdat);
            Ydat=1-Ydat;
        end
        if isfield(opts,'offset') && ~isempty(opts.offset)
            df=opts.offset*abs(Xdat(end)-Xdat(1));
            Xdat=Xdat+df;
        end
        % UNIQUE COMPARISON VALUES
        plot(Xdat,Ydat,sym,'color',color, 'markerface',markerface,'markersize',markersize,'LineWidth',LineWidth);
    end
    function xylim(obj,Xdat,opts)
        if nargin >= 3 || isempty(opts)
            if isfield(opts,'bXLim')
                bXLim=opts.bXLim;
            elseif ~isfield(opts,'bXLim') || isempty(opts.bXLim)
                bXLim=true;
            end
        end
        if bXLim
            [Xdat,~]=obj.get_dataXY();
            d=max(Xdat) - min(Xdat);
            m=d*0.05;
            xlim([min(Xdat)-m max(Xdat)+m]);
        end
        ylim([0 1]);
    end
    function ylabel(obj)
        ylabel('Prop. Cmp Chosen');
        yt=num2cell(yticks);
        yt=cellfun(@(x) sprintf('%.1f',x),yt,'UniformOutput',false);
        yticklabels(yt);;
    end
    function rylabel(obj,mult,opts)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        if ~isfield(opts,'bText')
            opts.bText=true;
        end
        if ~opts.bText
            return
        end

        text=obj.statsText(mult,newline);

        yyaxis right;
        ylabel(text);
        yticks('');
        yticklabels('');
        set(gca,'ycolor','k');
        set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','left') ;
        yyaxis left;

    end
    function strs=stdstr(obj,mult, names, bUnit)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        if nargin < 3
            names=[];
        end
        if nargin < 4
            bUnit=true;
        end
        n=length(obj.stds);
        strs=cell(n,1);
        vals=Vec.row(mult).*obj.stds;
        for i = 1:n
            if isempty(names)
                if n==1
                    I='X';
                else
                    I=['X' num2str(i)];
                end
            elseif length(names)==n
                I=names{i};
            elseif iscell(names)
                I=names{i};
            else
                I=names;
            end
            if bUnit
                I=[I ' '];
            else
                I='';
            end
            strs{i}=[I Num.toStr(vals(i))];
        end
    end
    function cmplabel(obj,meas,units)
        if nargin < 1 || isempty(meas)
            meas='x';
        end
        if nargin < 2 || isempty(units)
            units='a.u.';
        end
        meas(1)=upper(meas(1));
        units=['(' units ')'];
        xlabel([ meas ' ' units ]);
    end
    function cmpticks(obj,mult,xfrmt,Xdat,lbls)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        if nargin < 3
            xfrmt=[];
        end
        if nargin < 4 || isempty(Xdat)
            [Xdat,~]=obj.get_dataXY();
        end
        if nargin < 5 || isempty(lbls)
            lbls=xticklabels;
        elseif isnumeric(lbls)
            lbls=cellfun(@num2str,num2cell(lbls),'UniformOutput',false);
        end
        xticks(Xdat);
        if mult==1 && isempty(xfrmt)
            return
        end
        if isempty(xfrmt)
            spl=strsplit(lbls{1},'.');
            n=numel(spl{2});
            xfrmt=['%.' num2str(n) 'f'];
        end
        if ~isempty(lbls)
            lbls=cellfun(@(x) num2str(str2double(x)*mult,xfrmt),lbls,'UniformOutput',false);
            xticklabels(lbls);
        end
    end
    function obj = plot_boot_curve(obj,opts)

        if ~isfield(opts,'bExtendCurve')
            opts.bExtendCurve=obj.bExtendCurve;
        end
        if ~isfield(opts,'bAbs')
            bAbs=false;
        else
            bAbs=opts.bAbs;
        end

        a=min(obj.cmpX);
        b=max(obj.cmpX);
        if opts.bExtendCurve
            rang=b-a;
            a=a-rang*opts.bExtendCurve;
            b=b+rang*opts.bExtendCurve;
        end
        if isempty(obj.nX)
            obj.nX=200;
        end

        Xfit = linspace(a,b,obj.nX);

        % PSYCHOMETRIC FUNCTION FIT MEAN
        Yfit = obj.genGauss(Xfit,obj.muFit,obj.sigmaFit,obj.betFit,obj.DPcrt,obj.nIntrvl);
        if  obj.bExtendCurve
            ind=Yfit < .999 & Yfit > .001;
            Yfit=Yfit(ind);
            Xfit=Xfit(ind);
        end

        obj.get_bP();

        color=[];
        alpha=[];
        LineWidth=[];
        LineStyle={};
        if nargin >=2 && ~isempty(opts)
            if isfield(opts,'FillAlpha')
                alpha=opts.FillAlpha;
            end
            if isfield(opts,'LineWidth')
                LineWidth=opts.LineWidth;
            end
            if obj.bPlotCI && isfield(opts,'FillColor')
                color=opts.FillColor;
            elseif isfield(opts,'LineColor')
                color=opts.LineColor;
            end
            if isfield(opts,'LineStyle')
                LineStyle={opts.LineStyle};
            end
        else
        end
        if isempty(alpha)
            alpha=obj.CIalpha;
        end
        if isempty(LineWidth)
            LineWidth=obj.LineWidth;
        end
        if isempty(color) && obj.bPlotCI
            color=obj.color;
        elseif isempty(color)
            color=obj.lcolor;
        end
        if isempty(LineStyle)
            LineStyle={};
        end

        if obj.bPlotCI
            if obj.bP(1)
                mlh(1)=obj.muCI(1);
                mlh(2)=obj.muCI(2);
            else
                mlh(1)=obj.muFix;
                mlh(2)=obj.muFix;
            end
            if obj.bP(2)
                slh(1)=obj.sigmaCI(1);
                slh(2)=obj.sigmaCI(2);
            else
                slh(1)=obj.sigmaFix;
                slh(2)=obj.sigmaFix;
            end
            if obj.bP(2)
                blh(1)=obj.betCI(1);
                blh(2)=obj.betCI(2);
            else
                blh(1)=obj.betFix;
                blh(2)=obj.betFix;
            end

            c=Set.distribute(mlh,slh,blh);
            c=unique(c,'rows');
            Y=zeros(size(c,2),length(Xfit));
            for i = 1:size(c,1)
                ind=c(i,:);
                % m s b
                Y(i,:)  = obj.genGauss(Xfit,ind(1),ind(2),ind(3),obj.DPcrt,obj.nIntrvl);
            end
            Ymin=min(Y,[],1);
            Ymax=max(Y,[],1);
            obj.Ymin=min(Yfit);
            obj.Ymax=max(Yfit);
            obj.Xmin=min(Xfit);
            obj.Xmax=max(Xfit);

            %plot(Xfit,Ymin,'r'); hold on
            %plot(Xfit,Ymax,'r');
            %HERE
            patch([Xfit fliplr(Xfit)], [Ymin fliplr(Ymax)], color,'FaceAlpha',alpha,'EdgeColor','none'); hold on;
        end
        if bAbs && all(Xfit < 0)
            Xfit=abs(Xfit);
            Yfit=1-Yfit;
        end
        if isfield(opts,'offset') && ~isempty(opts.offset)
            [Xdat,Ydat]=obj.get_dataXY();
            df=opts.offset*abs(Xdat(end)-Xdat(1));
            Xfit=Xfit+df;
        end

        hold on;
        plot(Xfit,Yfit,LineStyle{:},'color',color,'LineWidth',LineWidth); hold on;


        % Data

    end
    function obj=format_boot_curve(obj,titl,units)
        if nargin < 2 || isempty(titl)
            titl='';
        else
            titl=[titl newline];
        end
        if nargin < 4 || isempty(units)
            units=obj.Xname;
        end
        if strcmp(units,'none')
            units='';
        end

        if obj.bTitle
        else
            titl=[];
        end
        %Axis.format(units,'Proportion Cmp Chosen',titl);
        Axis.format(units,'Prop. Cmp Chosen',titl);
        axis square;
    end

    function out=statsText(obj,mult,sep)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        if nargin < 3 || isempty(sep)
            %sep=newline;
            sep=', ';
        end
        out=[ ...
            'n=' num2str(size(obj.stdX,1)) sep ...
            '\mu='    num2str(mult*obj.muFit,'%2.2f') sep...
            '\sigma=' num2str(mult*obj.sigmaFit,'%2.2f') sep...
            'T=' num2str(mult*obj.tFit,'%2.2f') sep...
            '\beta='  num2str(obj.betFit,'%2.2f') sep ...
        ];
    end
    function obj = plot_boot_DP(obj)
        cmpXunq = transpose(unique(obj.cmpX));
        mm=[min(cmpXunq) max(cmpXunq)];

        c=polyfit(transpose(cmpXunq),obj.DPfit,1);
        f=@(x) c(1)*x + c(2);
        x=linspace(mm(1),mm(2),100);
        y=f(x);
        plot(x,y,obj.lcolor,'LineWidth',obj.LineWidth); hold on;

        plot(cmpXunq,obj.DPfit,[obj.shape obj.lcolor],'markersize',obj.markersize,'LineWidth',obj.LineWidth,'markerface',obj.markerface);

        titl=[ 'T='        num2str(obj.tFit,'%2.2f') ', N=' num2str(numel(obj.RCmpChs)) ];
        Axis.format(obj.Xname,'d''',titl);
        axis square;
        hold off;
    end
    function [] = plot_boot_T(obj)
        plot(obj.stdX,obj.tMU,[obj.shape obj.lcolor],'LineWidth',obj.LineWidth); hold on;
        hold off;
    end
    function [] = errorbar_boot_T(obj)
        %? tCI need to be subtracted?
        errorbar(obj.stdX,obj.tMU,obj.tCI,[obj.shape obj.lcolor],'LineWidth',obj.LineWidth); hold on;
        hold off;
    end
    function [] = plot_boot_mu_p(obj)
        obj.format_param_plot('\mu',obj.muSM,obj.muCI,obj.muFitDstb);
    end
    function [] = plot_boot_sigma_p(obj)
        obj.format_param_plot('\sigma',obj.sigmaSM,obj.sigmaCI,obj.sigmaFitDstb);
    end
    function [] = plot_boot_bet_p(obj)
        obj.format_param_plot('\beta',obj.betSM,obj.betCI,obj.betFitDstb);
    end
    function [] = plot_boot_T_p(obj)
        obj.format_param_plot('T',obj.tMU,obj.tCI,obj.tFitDstb);
    end
    function [] = format_param_plot(obj,name,mu,CI,dst)
        [H,B] = hist(dst,21);
        plot(B,H,'color',obj.lcolor,'LineWidth',obj.LineWidth);
        Axis.format(name,'Num Samples',['m=' num2str(mu,  '%2.2f') ', ' num2str(obj.CIsz) '%=[' num2str(CI(1), '%2.2f') ',' num2str(CI(2),'%2.2f')  ']']);
        lim=[min(B) max(B)];
        l=max(abs(mu-lim));
        xlim([mu-l mu+l]);
        %Text(.1,.9,{[num2str(obj.prcntUse) '% Data Used']},'ratio',18);
        axis square;
        hold off;
    end
    function obj=get_bP(obj)
        obj.bP=[isempty(obj.muFix) isempty(obj.muFix)  isempty(obj.betFix)];
    end

    %% IND PLOT
    function [] = plot_ind_fit(obj,PCdta)
        % PLOT FIT (IN HI-RES)
        [RCmpChsSel,cmpXsel,stdXsel]=obj.get_data();
        XcmpPlt = linspace(min(cmpXsel),max(cmpXsel),obj.nX);
        [PCplt,T]=obj.genGauss(XcmpPlt,obj.muFit,obj.sigmaFit,obj.betFit,obj.nIntrvl); hold on;

        cmpXUnq = unique(obj.cmpX)';

        plot(XcmpPlt,PCplt,'color',obj.lcolor,'linewidth',1.5); hold on

        % RAW DATA- COMPUTE PERCENT COMPARISON CHOSEN
        [PCdta] = PsyCurve.getpc(obj.sel);

        plot(cmpXUnq,PCdta,obj.shape,'color',obj.lcolor,'Linewidth',obj.LineWidth,'markersize',obj.markersize,'markerface',obj.markerface);

        % WRITE STUFF TO SCREEN
        Axis.writeText(1-.1,.1,{['n=' num2str(numel(RCmpChsSel))]},'ratio',18,'right');
        Axis.format('','',['T=' num2str(T,'%.2f') ': \mu=' num2str(obj.muFit,'%1.2f') ',\sigma=' num2str(obj.sigmaFit,'%1.2f') ',\beta=' num2str(obj.betFit,'%1.2f')]);
        xlim([min(obj.cmpXsel) max(cmpXsel)]+[-.1 .1]); ylim([0 1]);
        axis square;
        hold off;
    end
    function [] = plot_genGauss(obj,T)
        % XXX
        [RCmpChsSel,cmpXsel,stdXsel]=obj.get_data();
        plot(cmpXSel,PC,'color',color,'linewidth',2); hold on;
        Axis.format('','',['T=' num2str(T,'%.2f') ': \mu=' num2str(obj.muFit,'%1.2f') ',\sigma=' num2str(obj.sigmaFit,'%1.2f') ',\beta=' num2str(obj.betFit,'%1.2f') ',nIntrvl=' num2str(obj.nIntrvl)]);
        xlim(minmax(cmpXsel)+[-.1 .1]);
        ylim([0 1]);
        axis square;

        % WRITE PARAMETER VALUES AND N SMP TO SCREEN
        Axis.writeText(.075,.900,{['d''= ( | x - \mu| / \sigma )^{\beta}' ]},'ratio',20);
        Axis.writeText(.075,.775,{['d''_{crit}=' num2str(obj.DPcrt,'%.2f')]},'ratio',20);

        color='k';

        % THRESHOLD LINES
        plot(obj.muFit*[1 1],        [0                  0.5],'color',obj.lcolor,'linewidth',1);
        plot((obj.muFit+T)*[1 1],    [0     normcdf(0.5.*sqrt(obj.nIntrvl).*obj.DPcrt)],'color',obj.lcolor,'linewidth',1);
        plot([min(xlim) obj.muFit],  [0.5                0.5],'color',obj.lcolor,'linewidth',1);
        plot([min(xlim) obj.muFit+T],[normcdf(0.5.*sqrt(obj.nIntrvl).*obj.DPcrt)*[1 1]],'color',obj.lcolor,'linewidth',1);
        hold off;
    end
end
end
