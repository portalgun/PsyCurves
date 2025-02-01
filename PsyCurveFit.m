classdef PsyCurveFit < handle
% TODO allow beta fit for each cmp
properties
    Data

    % fix
    muFix
    sigmaFix
    betFix
    mFix
    yFix

    % fitInd
    muFitInd
    sigmaFitInd
    betFitInd
    mFitInd
    yFitInd

    bMuStdFix

    bLinear
    bLogLinear
    bNegLinear

    DPCrt
    nIntrvl

    bBoot
    nBoot
    bBootEachCmp
    bBest
    nBest
    CIsz
    prcntUse
    bootSeed

    fminOpts
    minFuncType

    best
    % TODO RM
    muFit
    sigmaFit
    betFit
    mFit
    yFit
    tFit
    %
    PCFit
    DPFit
    negLL

    % BOOT
    boot
    % STANDARD MEAN
    % TODO SM
    muSM
    sigmaSM
    betSM
    tSM
    % STANDARD ERROR
    % TODO SE
    muSE
    sigmaSE
    betSE
    tSE
    %
    % TODO CI
    muCI
    sigmaCI
    betCI
    tCI

    nLLFun

    measure='X'
    units
end
properties(Hidden=true)
    FitI
    Fix
    Any
    FitIP
    BFit
    NFit

    sel


    titl
    bP
    i

    stdX
    stds
    cmpX
    RCmpChs

    flds

end
methods(Static,Hidden)
    function P=getP()
        P={...
            'muFix',[],'Num.is';
            'sigmaFix',[],'Num.is';
            'yFix',[],'Num.is';
            'mFix',[],'Num.is';
            'betFix',1,'Num.is';
            ...
            'muFitInd',[],'Num.is';
            'sigmaFitInd',[],'Num.is';
            'yFitInd',[],'Num.is';
            'mFitInd',[],'Num.is';
            'betFitInd',[],'Num.is';
            ...
            'bMuStdFix',[],'Num.is';
            'bLinear',[],'Num.is';
            'bLogLinear',0,'Num.is';
            'bNegLinear',0,'Num.is';
            ...
            'DPCrt',1,'Num.is'; % XXX 1.36?
            'nIntrvl',2,'Num.isInt';
            ...
            'nBoot',[],'Num.isInt';
            'bBoot',[],'Num.is';
            'CIsz',68,'Num.is';
            'prcntUse',100,'Num.is';
            'minFuncType','fmincon','ischar';
            ...
            'bBest',1,'isbinary';
            'bBootEachCmp',[],'isbinary';
            'nBest',[],'Num.is';
            'bootSeed',1,'Num.is';
            ...
            'units','','ischar';
            'measure','','ischar';
        };
    end
end
methods(Static)
    function out=getDefaults()
        P=PsyCurveFit.getP();
        out=P(:,1:2)';
        out=struct(out{:});
    end
    function obj=new(stdX,cmpX,RCmpChs,varargin)
        if nargin < 1
            stdX=[];
        end
        if nargin < 2
            cmpX=[];
        end
        if nargin < 3
            RCmpChs=[];
        end
        Data=PsyCurveData(stdX,cmpX,RCmpChs);
        obj=PsyCurveFit(Data,varargin{:});
    end
end
methods
    function obj=PsyCurveFit(psyCurveData,varargin)
        if nargin < 1
            return
        end
        obj.Data=psyCurveData;

        if nargin < 2 || isempty(varargin)
            Opts=struct;
        elseif length(varargin)==1 && isstruct(varargin{1})
            Opts=varargin{1};
        else
            Opts=struct(varargin{:});
        end
        obj.parse_opts(Opts);
        obj.get_default_fminOpts();
    end
    function get_default_fminOpts(obj)
        % SET FMINCON OPTIONS
        if strcmp(obj.minFuncType,'fmincon')
            obj.fminOpts             = optimset('fmincon');
            obj.fminOpts.MaxPCGIter=[];
            obj.fminOpts.Algorithm   = 'sqp';
            obj.fminOpts.LargeScale  = 'off';
            obj.fminOpts.UseParallel = 'never';
            obj.fminOpts.Display     = 'none';
            %obj.fminOpts.Display='iter';
            obj.fminOpts.MaxIter     = 500;
            obj.fminOpts.MaxFunEvals = 500;
            %obj.fminOpts.TolCon      = 500;
        elseif strcmp(obj.minFuncType,'fminsearch')
            obj.fminOpts             = optimset('fminsearch');
            obj.fminOpts.UseParallel = 'never';
            obj.fminOpts.Display     = 'off';
            obj.fminOpts.MaxIter     = 500;
        end
    end
    function obj=run(obj,bRefit)
        obj.init_params();

        if nargin < 2
            bRefit=false;
        end
        if bRefit
            lMFit=obj.muFit;
            lTFit=obj.tFit;
            lSFit=obj.sigmaFit;
            lBFit=obj.betFit;
            lNegLL=obj.negLL;
        end
        if obj.bBest
            obj.fit_best();
        else
            obj.fit_basic();
        end
        if bRefit && obj.negLL < lNegLL
            obj.muFit  =lMFit;
            obj.tFit  =lTFit;
            obj.sigmaFit  =lSFit;
            obj.betFit  =lBFit;
            obj.negLL =lNegLL;
        end

        if obj.bBoot
            obj.fit_boot();
        end
    end
%- Plot
    function plot_negLL(obj)
        [RCmpChs,cmpX,stdX]=obj.sel_data();
        m=unique(stdX);
        s=.001:.001:.02;
        b=1;

        % SET INITIAL PARAMETER VALUES
        nl=zeros(length(s),1);
        m0  = obj.muFix;
        s0  = obj.sigmaFix;
        b0  = obj.betFix;
        if isempty(m0); m0 = mean([min(cmpX) max(cmpX)]); m0 = m0  + .1.*randn; end
        if isempty(s0); s0 = diff(Num.minMax(abs(cmpX)))./6;        s0 = s0  + .1.*s0.*randn; end
        if isempty(b0); b0 = 1;                                 b0 = b0  + .1.*b0.*randn; end
        p0 = [m0 s0 b0];

        fun=@(p) PsyCurveFit.negLLFunc(p,cmpX,RCmpChs,obj.DPCrt,obj.nIntrvl,obj.muFix,obj.sigmaFix,obj.betFix);
        for i = 1:length(s)
            nl(i)=fun([p0(1),s(i),p0(3)]);
        end
        plot(s,nl);
    end
%- Text
    function summary(obj)
        dispV(obj.muFit,'mu')
        dispV(obj.sigmaFit,'sigma')
        dispV(obj.tFit,'t')
        dispV(obj.mFit,'m')
        dispV(obj.yFit,'y')
        dispV(obj.betFit,'bet')
        dispV(obj.negLL,'negLL')
    end
%- Plot
    function plot(obj,varargin)

        figure(1)
        hold off;

        X=obj.getX();
        Y=obj.gen_gauss_X(X)*100;

        n=size(Y,1);
        lcolors=zeros(n,3);
        for i = 1:n
            plot(X,Y(i,:),varargin{:},'LineWidth',2);
            lcolors(i,:)=lastColor();
            hold on;
        end

        %subPlot([1 2],1,2);
        obj.Data.plotPC('LineStyle','none','Colors',lcolors,'Marker','o',varargin{:});
        axis square;


        if ~isempty(obj.units)
            ulbl=[' (' obj.units ')'];
        else
            ulbl='';
        end
        xlabel([obj.measure ulbl]);
        ylabel('Percent Comparison Chosen');
        Axis.format();
    end
    function plotT(obj,varargin)
        if nargin > 1
            % TODO varargin
            TODO varargin
        end
        figure(2)
        hold off;

        % data
        if obj.bLinear
            lstyle='none';
        else
            lstyle='-';
        end

        Xd=obj.Data.stdXUnq;
        Yd=obj.tFit;

        plot(Xd,Yd,'ok','LineStyle',lstyle,'MarkerFaceColor','k');
        if ~isempty(obj.units)
            ulbl=[' (' obj.units ')'];
        else
            ulbl='';
        end
        ylabel([' Threshold ' ulbl]);
        xlabel([obj.measure ulbl]);
        Axis.format();
        if ~obj.bLinear
            return
        end

        % line
        X=obj.getX();
        [Ys,Yt]=PsyCurve.lin2sigma(X,obj.mFit,obj.yFit,obj.betFit(1),obj.bLogLinear,obj.DPCrt);

        hold on;
        plot(X,Yt,'k','LineStyle','-','MarkerFaceColor','none');
        if obj.bLogLinear
            set(gca,'YScale','log');
        end



    end
end
methods(Hidden)
    function obj=parse_opts(obj,Opts)
        if ~exist('Opts','var')
            Opts=struct;
        end

        Args.parse(obj,obj.getP,Opts);

        % Linear
        if isempty(obj.bLinear)
            obj.bLinear=obj.bLogLinear;
        elseif obj.bLogLinear && ~obj.bLinear
            error('bLogLinear and bLinear are conflicting')
        end

        % Best
        if isempty(obj.bBest)
            if isempty(obj.nBest)
                obj.bBest=false;
            else
                obj.bBest = obj.nBest > 1;
            end
        end
        if isempty(obj.nBest) && obj.bBest
            obj.nBest=50;
        end

        % Boot
        if isempty(obj.bBoot)
            if isempty(obj.nBoot)
                obj.bBoot=false;
            else
                obj.bBoot = obj.nBoot > 1;
            end
        end
        if isempty(obj.nBoot) && obj.bBoot
            obj.nBoot=1000;
        elseif isempty(obj.nBoot)
            obj.nBoot=0;
        end
        if isempty(obj.bBootEachCmp)
            obj.bBootEachCmp=~obj.bLinear;
        end

        % bMuStdFix
        if isempty(obj.bMuStdFix)
            obj.bMuStdFix=obj.bLinear && ~isempty(obj.muFix);
        end

    end
%- MAIN
    function fit_basic(obj)
    %NON BOOTSTRAPPED FIT
        obj.sel_data();
        S = obj.init_sel(1);
        S = obj.fit(S);

        obj.pack_basic(S);
    end
    function fit_best_linear(obj)
        % XXX
    end
    function fit_best(obj)
        obj.sel_data();
        S=obj.init_sel(obj.nBest);

        for i = 1:obj.nBest
            S(i) = obj.fit(S(i));
            S(i) = obj.gen_gauss_cmp(S(i));
        end
        obj.pack_best(S);
    end
    function fit_boot(obj)
        obj.sel_boot_data(0,[]);
        S=obj.init_sel(obj.nBoot);

        rng(obj.bootSeed,'twister');
        sds=randi(2^32-1,obj.nBoot,1);

        for i = 1:obj.nBoot
            %obj.i=i;
            sd=sds(i);

            % Sample
            obj.sel_boot_data(sd,[]);
            S(i) = obj.fit(S(i));
            S(i) = obj.gen_gauss_cmp(S(i));
        end
        obj.pack_boot(S);
    end
%- Init
    function init_params(obj)
        stdXU=unique(obj.Data.stdX);
        nStd=numel(stdXU);
        if obj.bMuStdFix
            obj.muFix=stdXU;
        end
        % fld, n, bUse, bUniform
        flds={...
            'mu',   nStd, 1,           0;
            'sigma',nStd,~obj.bLinear, 0;
            'y',    1,    obj.bLinear, 1;
            'm',    1,    obj.bLinear, 1;
            'bet',  nStd, 1,           1;
        };

        % XXX make sure sigma OR (y & m)
        count=0;
        Any=struct();
        FitIP=struct();
        Fix=struct();
        FitI=struct();
        BFit=struct();
        for i = 1:length(flds)
            fld=flds{i,1};
            n=flds{i,2};
            bUse=flds{i,3};
            bUniform=flds{i,4};

            if ~bUse
                continue
            end

            fInd=[fld 'FitInd'];
            fFix=[fld 'Fix'];

            FitI.(fld)=obj.(fInd);
            Fix.(fld)=obj.(fFix);

            if ~isempty(Fix.(fld))
                if numel(Fix.(fld))==1
                    Fix.(fld)=repmat(Fix.(fld),n,1);
                end
            else
                Fix.(fld)=nan(n,1);
            end

            % fitIndeces
            if all(isnan(Fix.(fld))) && isempty(FitI.(fld))
                %  nothing specified -> assumes fit
                if bUniform
                    FitI.(fld)=ones(n,1);
                else
                    FitI.(fld)=(1:n)';
                end
            elseif ~isempty(FitI.(fld)) && any(~isnan(Fix.(fld)) & FitI.(fld) > 0)
                % fix and fInd specified and conflicting
                error([ fFix ' and ' fInd ' are conflicting'])
            elseif ~isempty(FitI.(fld))
                %  fInd specified
                nIn=numel(FitI.(fld));
                if nIn == 1
                    FitI.(fld)=repmat(FitI.(fld),n,1);
                elseif nIn~=n
                    error(['Number of elements in field ' fInd ' must be ' num2str(n) '.']);
                elseif Vec.isRow(FitI.(fld))
                    FitI.(fld)=FitI.(fld)';
                end
            end

            % BFit
            BFit.(fld)=isnan(Fix.(fld));

            % Any
            Any.(fld)=any(BFit.(fld));

            % NFit
            NFit.(fld)=sum(unique(FitI.(fld))>0);

            % fitInd param
            if Any.(fld)
                FitIP.(fld)=count+(1:NFit.(fld));
                count=count+NFit.(fld);
            end
        end
        obj.flds=fieldnames(Any);

        obj.Fix=Fix;
        obj.FitI=FitI;
        obj.BFit=BFit;
        obj.Any=Any;
        obj.NFit=NFit;
        obj.FitIP=FitIP;

    end
    function varargout=get_params(obj,param)
        N=size(obj.flds,1);
        varargout=cell(1,N);
        for i = 1:N
            fld=obj.flds{i,1};
            varargout{i}=obj.Fix.(fld);
            if obj.Any.(fld)
                p=param(obj.FitIP.(fld));
                varargout{i}(obj.BFit.(fld))=p(obj.FitI.(fld));
            end
        end
    end
    function S=init_sel(obj,n)

        nC=0;
        nS=obj.sel.nS;
        cS=obj.sel.cS;
        for is = 1:nS
            ind=cS==is;
            nc=numel(unique(obj.sel.cmpX(ind)));
            if nc > nC
                nC=nc;
            end
        end
        z=nan(nS,1);
        zc=nan(nS,nC);

        if ~isfield(obj.flds,'sigma')
            ym=nan;
        else
            ym=[];
        end

        S=struct(...
                 'negLL',nan,...
                 ...
                 'mu',z,...
                 'sigma',z,...
                 'y',ym,...
                 'm',ym,...
                 'bet',nan,...
                 't',z,...
                 ...
                 'PCFit',zc,...
                 'DPFit',zc);
        S=repmat(S,n,1);
    end
%- SEL DATA
    function sel_data(obj)
        obj.sel=obj.Data;
    end
    function sel_boot_data(obj, sd,nSmp)
        obj.sel=obj.Data.selBoot(obj.bBootEachCmp,sd,nSmp,obj.prcntUse);
    end
%- PACK
    function pack_basic()
        obj.pack_dist(S);
        obj.select_from_dist(1);
        obj.muFit=vertcat(S.muFit);
        obj.sigmaFit=vertcat(S.sigmaFit);
        obj.betFit=vertcat(S.betFit);
        obj.tFit=vertcat(S.tFit);
        obj.negLL=vertcat(S.negLL);
    end
    function pack_best(obj,S)
        obj.pack_dist(S,'best');
        ind=find(obj.best.negLL==min(obj.best.negLL),1,'first');
        obj.select_from_dist(ind,'best');
    end
    function obj=pack_boot(obj,S)
        obj.pack_dist(S,'boot');


        CIlohi = 0.5*(1-obj.CIsz/100) + [0 obj.CIsz/100];

        flds={'mu','sigma','bet','t'};
        for i = 1:length(flds)
            fld=flds{i};
            dat=obj.boot.(fld);
            obj.([fld 'CI'])=quantile(dat,CIlohi,3);
            obj.([fld 'SE'])=    std(dat,[],3);
            obj.([fld 'SM'])=   mean(dat,3);
        end

        return

        % CONFIDENCE INTERVAL: LO & HI BOUNDS
        % BOOSTRAPPED CONFIDENCE INTERVALS
        obj.muCI    = quantile(fld.muFitDstb(~isnan(fld.muFitDstb)), CIlohi);
        %obj.sigmaCI = quantile(fld.sigmaFitDstb(~isnan(fld.sigmaFitDstb)), CIlohi);
        obj.betCI   = quantile(fld.betFitDstb(~isnan(fld.betFitDstb)), CIlohi);
        obj.tCI     = quantile(fld.tFitDstb(~isnan(fld.tFitDstb)), CIlohi);
        % BOOSTRAPPED STD ERR OF STATISTIC
        obj.muSE    = std(fld.muFitDstb(~isnan(fld.muFitDstb)));
        obj.sigmaSE = std(fld.sigmaFitDstb(~isnan(fld.sigmaFitDstb)));
        obj.betSE   = std(fld.betFitDstb(~isnan(fld.betFitDstb)));
        obj.tSE     = std(fld.tFitDstb(~isnan(fld.tFitDstb)));
        % BOOSTRAPPED MEAN
        obj.muSM    = mean(fld.muFitDstb(~isnan(fld.muFitDstb)));
        obj.sigmaSM = mean(fld.sigmaFitDstb(~isnan(fld.sigmaFitDstb)));
        obj.betSM   = mean(fld.betFitDstb(~isnan(fld.betFitDstb)));
        obj.tMU     = mean(fld.tFitDstb(~isnan(fld.tFitDstb)));
    end
    function obj=pack_dist(obj,S,name)
        flds=fieldnames(S);
        for i = 1:length(flds)
            fld=flds{i};
            %if any(strcmp(fld,{'PCFit','DPFit'}))
            %    % XXX
            %    continue
            %end
            %S(1).sigma

            obj.(name).(fld)=cat(3,S.(fld));
        end

    end
%- SELECT OUTPUT
    function select_from_dist(obj,ind,name)
        flds=fieldnames(obj.(name));
        for i = 1:length(flds)
            fld=flds{i};
            if any(strcmp(fld,{'negLL','PCFit','DPFit'}))
                ofld=fld;
            else
                ofld=[fld 'Fit'];
            end
            obj.(ofld)=obj.(name).(fld)(:,:,ind);
        end
    end
    function select_mean(obj)
        % XXX
        obj.muFit = obj.muSM;
        obj.sigmaFit = obj.sigmaSM;
        obj.betFit = obj.betSM;

        [RCmpChs,cmpX,stdX]=obj.sel_data();
        p0 = [obj.muFit obj.sigmaFit obj.betFit];
        obj.negLL=negLLFunc(p,cmpX,RCmpChs,obj.DPCrt,obj.nIntrvl);

        obj.tFit = obj.tMU;
    end
    function select_err(obj)
        % XXX
        [Xdat,Ydat]=obj.sel_dataXY();
        Ydat=Vec.row(Ydat);
        ctr=ceil(length(Xdat/2));
        %Ydat(ctr)=[];
        %Xdat(ctr)=[];

        E2=zeros(obj.nBoot,1);
        for i = 1:length(obj.nBoot)
            Y = PsyCurve.genGauss(Xdat,obj.muFitDstb(i,:),obj.sigmaFitDstb(i,:),obj.betFitDstb(i,:));
            E2(i)=mean((Y-Ydat).^2);
        end
        ind=find(E2==min(E2),1,'first');

        obj.select_from_dist(ind);
    end
%- FIT
    function s = fit(obj,s)

        obj.getNegLLFun();
        [pLB,pUB,p0]=obj.getBounds();

        if obj.bLinear
            A=[];
            b=[];
            Ai=[];
            bi=[];
            if obj.bLogLinear
                nlfun=@(p) obj.nlconstr(p);
            else
                nlfun=[];
            end
        else
            A=[];
            b=[];
            Ai=[];
            bi=[];
            nlfun=[];
        end

        % MINIMIZE NEGATIVE LOG-LIKELIHOOD
        switch obj.minFuncType
        case 'fmincon'
            [pFit,s.negLL] = fmincon(obj.nLLFun,p0,Ai,bi,A,b,pLB,pUB,nlfun,obj.fminOpts);
        case 'fminsearch'
            [pFit,s.negLL] = fminsearch(neg_LL,p0,obj.fminOpts);
        otherwise
            error()
        end

        % FINAL FIT PARAMETERS
        n=length(obj.flds);
        args=cell(1,n);
        [args{:}]=obj.get_params(pFit);
        for j = 1:n
            fld=obj.flds{j};
            s.(fld)=args{j};
        end

        % Handle equivalent solutions
        if ~ismember('sigma',obj.flds)
            [s.sigma,s.t]=PsyCurve.lin2sigma(s.mu,s.m,s.y,s.bet,obj.bLogLinear,obj.DPCrt);
            if any(sign(s.sigma) ~=  sign(s.t))
                % - -
                s.m=-s.m;
                s.y=-s.y;
                [~,s.t]=PsyCurve.lin2sigma(s.mu,s.m,s.y,s.bet,obj.bLogLinear,obj.DPCrt);
                if any(sign(s.sigma) ~=  sign(s.t))

                    % - +
                    s.m=-s.m;
                    [~,s.t]=PsyCurve.lin2sigma(s.mu,s.m,s.y,s.bet,obj.bLogLinear,obj.DPCrt);

                    if any(sign(s.sigma) ~=  sign(s.t))

                        % + -
                        s.m=-s.m;
                        s.y=-s.y;
                        [~,s.t]=PsyCurve.lin2sigma(s.mu,s.m,s.y,s.bet,obj.bLogLinear,obj.DPCrt);
                    end
                end
            end
        end

    end
    function [C,Ceq]=nlconstr(obj,p)
        Ceq=[];

        [mu,y,m,bet]=obj.get_params(p);
        X=Num.minMax([obj.Data.cmpXUnq(:);mu])';
        [sigma,t]=PsyCurve.lin2sigma(X,m,y,bet(1),obj.bLogLinear,obj.DPCrt);
        %t=PsyCurve.sigma2thresh(sigma,obj.DPCrt,bet(1));
        tol=10^-6;
        %if any(isnan(t))
        %    m
        %    y
        %    sigma
        %end
        %C=[-sigma+tol; -t+tol; t-50; -(y*m)+tol];
        C=[-t+tol];
    end
%- Bounds & P0
    function [pLB,pUB,p0]=getBounds(obj)
        cmpXUnq=obj.sel.cmpXUnq;

        d=2*(cmpXUnq-mean(cmpXUnq,2))+mean(cmpXUnq,2);
        mu0=[min(d,[],2) max(d,[],2)];

        mm=[max(cmpXUnq,[],2) min(cmpXUnq,[],2)];
        d=abs(diff(mm,[],2));
        sigma0=[0.02 2.00].*d;

        mu=[];
        sigma=[];
        bet=[];
        m=[];
        y=[];
        if isfield(obj.Any,'mu') && obj.Any.mu
            mu=mu0;
            mu(~obj.FitI.mu)=[];
        end
        if isfield(obj.Any,'sigma') && obj.Any.sigma
            sigma=sigma0;
            sigma(~obj.FitI.sigma)=[];
        end
        if isfield(obj.Any,'y') && obj.Any.y
            % TODO
            y=[-10 10];
            %if obj.bLogLinear
            %end
        end
        if isfield(obj.Any,'m') && obj.Any.m
            t=PsyCurve.sigma2thresh(sigma0,obj.DPCrt);
            t=max(t(:))/min(t(:));
            if obj.bLogLinear
                t=log(t);
            end
            m=[-t t];
        end
        if obj.Any.bet
            bet =  repmat([0.35 3.00],obj.NFit.bet,1);
        end
        p=cat(1,mu,sigma,y,m,bet);

        % Param0
        r=abs(diff(p,[],2));
        out=nan;
        c=0;
        best=inf;
        while true
            c=c+1;
            if c > 100
                error('Parameters do not appear to allow for valid evaluation')
            end
            p00=rand([size(p,1),1]).*r+p(:,1);
            out=obj.nLLFun(p00);
            if out < best
                best=out;
                p0=p00;
            end
            if c > 20 && (~isnan(best) && ~isinf(best))
                break
            end

        end

        pLB=p(:,1);
        pUB=p(:,2);

    end

%- NEGLL
    function getNegLLFun(obj)
        if obj.bLinear
            obj.nLLFun=@(p) obj.negLLFuncLinear_(p);
        else
            obj.nLLFun=@(p) obj.negLLFunc_(p);
        end
    end
    function out=negLLFunc_(obj,p)
        [mu,sigma,bet]=obj.get_params(p);
        if numel(mu)==1
            out=PsyCurveFit.negLLFunc(mu,sigma,bet,obj.sel.cmpX,obj.sel.RCmpChs, obj.nIntrvl);
        else
            out=PsyCurveFit.negLLFuncMult(mu,sigma,[],bet, obj.sel.stdX,obj.sel.cmpX,obj.sel.RCmpChs, obj.DPCrt,obj.nIntrvl,obj.bLogLinear,false);

        end
    end
    function out=negLLFuncLinear_(obj,p)
        [mu,y,m,bet]=obj.get_params(p);
        out=PsyCurveFit.negLLFuncMult(mu,y,m,bet, obj.sel.stdX,obj.sel.cmpX,obj.sel.RCmpChs, obj.DPCrt,obj.nIntrvl,obj.bLogLinear,true);
    end
    function Y=gen_gauss_X(obj,X);
        nS=obj.Data.nS;
        Y=zeros(nS,numel(X));
        for i=1:nS
            Y(i,:)=PsyCurve.genGauss(X,obj.muFit(i),obj.sigmaFit(i),obj.betFit(i),obj.nIntrvl);
        end
    end
    function s=gen_gauss_cmp(obj,s)
        if isempty(s.sigma)
            sigma=PsyCurve.lin2sigma(s.mu,s.m,s.y,obj.bLogLinear);
            %sigma=repmat(sigma,numel(s.mu),1);
        else
            sigma=s.sigma;
        end
        [stdXU,~,C]=unique(obj.sel.stdX);
        nC=numel(stdXU);
        z=zeros(nC,1);

        for ic = 1:nC
            ind=C==ic;
            cmpX=unique(obj.sel.cmpX(ind));

            [s.PCFit(ic,:),s.DPFit(ic,:)] = PsyCurve.genGauss(cmpX,s.mu(ic),sigma(ic),s.bet(ic),obj.nIntrvl);
        end
        if ~obj.bLinear
            s.t=PsyCurve.sigma2thresh(sigma,obj.DPCrt);
        end
    end
    function X=getX(obj,n)
        if nargin < 2 || isempty(n)
            n=1000;
        end
        rng=Num.minMax(obj.Data.cmpXUnq);
        X=linspace(rng(1),rng(2),n);
    end
end
%- NLL Funcs
methods(Static)
    function out=negLLFunc(mu,sigma,bet, cmpX,RCmpChs, nIntrvl);
        out=-(sum(log(     PsyCurve.genGauss(cmpX(RCmpChs==1),mu,sigma,bet, nIntrvl) )) + ...
              sum(log( 1 - PsyCurve.genGauss(cmpX(RCmpChs==0),mu,sigma,bet, nIntrvl) )) );
    end
    function out=negLLFuncMult(mu,y,m,bet, stdXAll,cmpXAll,RCmpChsAll, DPCrt,nIntrvl,bLog,bLinear)
        % mu, y, bet have numel = nStdXU

        if bLinear
            [sigma,t]=PsyCurve.lin2sigma(mu,m,y,bet,bLog,DPCrt);
        else
            sigma=y;
        end

        [stdXU,~,C]=unique(stdXAll);
        nC=numel(stdXU);
        nLL=zeros(nC,1);
        for ic = 1:nC
            ind=C==ic;
            nLL(ic)=PsyCurveFit.negLLFunc(mu(ic),sigma(ic),bet(ic), cmpXAll(ind),RCmpChsAll(ind), nIntrvl);
        end
        out=sum(nLL);
        if any(isinf(out)) || any(isnan(out))
            %out=inf;
            %out=nan;
        end
        %[m y]
        %t
        %sigma
        if any(isnan(out))
            dk
        end
    end
end
methods(Static,Hidden)
    function test0()
        nTrlPerCmp=1000;  % number of trials per comparison
        stdX=[7.5]; % standard
        cmpX=[5.1540    6.3300    7.5000    8.6700    9.8460]; % comparisons

        sigma=[0.8];
        cr=5;
        mu=7.5;

        %S=PsyCurveData.genData(nTrlPerCmp,stdX,cmpX,[],sigma);
        S=PsyCurveData.genData(nTrlPerCmp,stdX,cmpX,[],sigma,cr);
        Fig.new('dv')
        hist(S.DV)

        Fig.new('rv')
        plot(S.nRCmpChs)

        %% Set Fitting Options
        Opts=PsyCurveFit.getDefaults();
        Opts.nIntrvl=2;
        Opts.units='arcmin';
        Opts.measure='Disparity';
        %Opts.DPCrt=1;
        %Opts.nBest=5;
        %Opts.nBoot=10;
        Opts.muFix=stdX;


        %% Construct
        pc=PsyCurveFit.new(S.stdX,S.cmpX,S.RCmpChs,Opts);

        %% Run
        pc.run()
        Fig.new('plot')
        pc.plot()
    end
    function pc=test2()
        stdX=[3.75;     5.625;    7.5000;   9.375;     11.25]; % standards
        cmpX=[-2.3460 -1.1700          0   1.1700     2.3460];
        mu=[];
        sigma=[0.5 0.6 0.8 0.9 1.0]; % standard deviation
        nTrlPerCmp=100;  % number of trials per comparison

        S=PsyCurveData.genData(nTrlPerCmp,stdX,cmpX,mu,sigma);

        %% Set Fitting Options
        Opts=PsyCurveFit.getDefaults();
        Opts.nIntrvl=2;
        Opts.DPCrt=1;
        Opts.nBest=5;
        Opts.nBoot=10;
        Opts.muFix=stdX;
        Opts.bLinear=true;
        Opts.bBootEachCmp=true;
        Opts.measure='Disparity';
        Opts.units='arcmin';

        %% Construct
        pc=PsyCurveFit.new(S.stdX,S.cmpX,S.RCmpChs,Opts);

        %% Run
        pc.run();

        %% Plot Curves
        pc.plot();

        %% Plot Threshold
        pc.plotT();

    end
    function test1()
        % TODO
        %    test log linear
        %

        n=1;

        mu=[-2 -1 0 1 2 ];
        cr=0;
        nTrlPerCmp=100;

        %sigma=[1];
        sigma=[.5; 1; 1.5];

        %stdX=7.5;
        stdX=[3.75 7.5 11.25];
        %stdX=[11.25 7.5 3.75];

        Opts=PsyCurveFit.getDefaults();
        Opts.bLinear=false;
        Opts.DPCrt=1;
        Opts.nBest=1;

        Opts.bMuStdFix=true;
        Opts.bLinear=true;
        Opts.bLogLinear=true;

        %[S,G]=PsyCurveData.genNorm(mu,sigma,cr,nTrlPerCmp);
        %Dat=PsyCurveData(S.sC,S.cC,S.RCmpChs);

        [S,G]=PsyCurveData.gen(nTrlPerCmp,stdX,mu,sigma,cr);
        S
        dk
        Dat=PsyCurveData(S.stdX,S.cmpX,S.RCmpChs);
        Dat.G=G;

        obj=PsyCurveFit(Dat,Opts);
        obj.run();

        obj.plot();
        obj.plotT();
        obj.summary();


        assignin('base','obj',obj);

        %for i = 1:dat.nCmp
        %    plot(
        %end



    end
end
end

