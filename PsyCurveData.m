classdef PsyCurveData < handle
properties
    RCmpChs

    stdX
    stdXUnq
    cS
    nS

    cmpX
    cmpXUnq
    cC
    nC

    %lvlX
    %lvlXUnq
    cL
    nL

    n
    PC

    G
end
methods(Static)
end
methods
    function obj=PsyCurveData(stdX,cmpX,RCmpChs)
        if nargin < 1
            return
        end
        if nargin < 4
            bFlipCmpChs=[];
        end
        [stdX,cmpX,RCmpChs]=PsyCurveData.parse(stdX,cmpX,RCmpChs,bFlipCmpChs);
        obj.parse_lvls_(stdX,cmpX,RCmpChs);
        obj.getPC_();
    end
end
methods(Hidden)
    function parse_lvls_(obj,stdX,cmpX,RCmpChs)
        obj.stdX=stdX;
        obj.cmpX=cmpX;
        obj.RCmpChs=RCmpChs;

        obj.n=size(obj.cmpX,1);

        [obj.stdXUnq,~,obj.cS]=unique(stdX);
        obj.nS=numel(obj.stdXUnq);

        % CmpXUnq
        obj.cmpXUnq=cell(obj.nS,1);
        for is = 1:obj.nS
            indS=obj.cS==is;
            [obj.cmpXUnq{is}]=unique(obj.cmpX(indS))';
        end
        obj.cmpXUnq=vertcat(obj.cmpXUnq{:});
        obj.nC=size(obj.cmpXUnq,2);

        % cC & c
        obj.cC=zeros(obj.n,1);
        obj.cL=zeros(obj.n,1);
        c=0;
        for is = 1:obj.nS
            indS=obj.cS==is;
            for ic = 1:obj.nC
                c=c+1;
                indC=obj.cmpX==obj.cmpXUnq(is,ic);
                obj.cC(indC & indS)=ic;
                obj.cL(indC & indS)=c;
            end
        end
        obj.nL=c;
    end
    function getPC_(obj)
        obj.PC=zeros(obj.nS,obj.nC);
        for is=1:obj.nS
            indS=obj.cS==is;
            for ic=1:obj.nC
                indC=obj.cC==ic;
                obj.PC(is,ic)=sum(obj.RCmpChs(indC & indS)==1)/sum(indC & indS)*100;
            end
        end
    end
end
methods
    function out=copy(obj)
        out=PsyCurveData();
        out.init_(obj.stdX,obj.cmpX,obj.RCmpChs);
    end
    function struct(obj)
        out=struct('stdX',obj.stdX,'cmpX',obj.cmpX,'RCmpChs',obj.RCmpChs);
    end
    function out=selBoot(obj, bBootEachCmp, sd, nSmp, prcntUse)
        if nargin < 4
            nSmp=[];
        end
        if nargin < 5
            prcntUse=[];
        end
        if isempty(nSmp) && isempty(prcntUse)
            prcntUse=100;
        end

        S0=obj.struct();
        S=struct();
        if bBootEachCmp
            S.RCmpChs=cell(obj.nS,1);
            S.stdX=cell(obj.nS,1);
            S.cmpX=cell(obj.nS,1);
            for ic = 1:obj.cS

                ind=find(obj.cS==ic);
                if isempty(nSmp)
                    nSmp=round(numel(ind).*prcntUse./100);
                end
                rng('twister',sd);
                indSmp = ind(randi(nSmp,nSmp,1));

                S.RCmpChs{ic}=   S0.RCmpChs(indSmp);
                S.stdX{ic}   =   S0.stdX(indSmp);
                S.cmpX{ic}   =   S0.cmpX(indSmp);
            end
            S.RCmpChs=vertcat(S.RCmpChs{:});
            S.stdX=vertcat(S.stdX{:});
            S.cmpX=vertcat(S.cmpX{:});
        else
            if isempty(nSmp)
                nSmp=round(numel(C).*prcntUse./100);
            end
            indSmp = ind(randi(nSmp,nSmp,1));
            S.RCmpChs=S0.RCmpChs(indSmp);
            S.stdX=S0.stdX(indSmp);
            S.cmpX=S0.cmpX(indSmp);

        end
        out=PsyCurveData(S.stdX,S.cmpX,S.RCmpChs);
    end
    function plotPC(obj,varargin)
        [colors,varargin,bSuccess]=Args.getPair('Colors',varargin{:});
        for is = 1:obj.nS
            if bSuccess
                cargs={'MarkerEdgeColor',colors(is,:),'MarkerFaceColor',colors(is,:)};
            else
                cargs={};
            end
            plot(obj.cmpXUnq(is,:),obj.PC(is,:),cargs{:},varargin{:});
            hold on;
        end
        ylim([0 100]);
    end
end
methods(Static)
    function [stdX,cmpX,RCmpChs]=parse(stdX,cmpX,RCmpChs, bFlipCmpChs)
        if nargin < 4 || isempty(bFlipCmpChs)
            bFlipCmpChs=false;
        end
        % XXX make sure std cmp Rchs are same size
        [n,m]=size(stdX);
        if n < m && n > 1 && m > 1
            stdX=stdX';
        end
        stdX=rmUniformCols(stdX);
        cmpX=rmUniformCols(cmpX);
        if numel(stdX) == 1
            stdX = stdX*ones(numel(cmpX),1);
        end

        % sizes
        if size(stdX,2)    ~= 1
            stdX    = stdX(:);
        end
        if size(cmpX,2)    ~= 1
            cmpX    = cmpX(:);
        end
        if size(RCmpChs,2) ~= 1
            RCmpChs = RCmpChs(:);
        end

        % rm nan
        nnanind=~isnan(stdX) & ~isnan(cmpX) & ~isnan(RCmpChs);
        RCmpChs=logical(RCmpChs(nnanind));
        cmpX=cmpX(nnanind);
        stdX=stdX(nnanind);

        % vectorize
        stdX=stdX(:);
        cmpX=cmpX(:);
        RCmpChs=RCmpChs(:);

        % sort
        [~,ind]=sortrows([stdX cmpX]);
        stdX=stdX(ind);
        cmpX=cmpX(ind);
        RCmpChs=RCmpChs(ind);

        function X=rmUniformCols(X)
            [n,m]=size(X);
            if n==1
                return
            end
            ind=false(1,m);
            for i = 1:m
                ind(i)=all(X(:,i)==X(1,i));
            end
            if all(ind)
                ind(1)=false;
            end
            X(:,ind)=[];
        end
    end
    function out=gen2(mu,sigma,cr,nTrlPerCmp)
        nStd=numel(stdX);
        nCmp=size(cmpX,2);

        % mu
        if nargin < 1 || isempty(cmpX)
            mu=0;
        elseif size(mu,2)==1
            mu=repmat(mu,1,nMu);
        end

        % sigma
        if nargin < 2 || isempty(sigma)
            sigma=ones(1,nMu);
        elseif size(sigma,2)==1
            sigma=repmat(sigma,1,nMu);
        end

        % cr
        if nargin < 3
            cr=zeros(1,nMu);
        elseif size(cr,2)==1
            cr=repmat(cr,1,nMu);
        end
        if nargin < 4
            nTrlPerCmp=100;
        end

        S=struct();
        for i =1:nStd
            S(i)=genCmps(stdX(i),mu(i,:),sigma(i,:),cr(i,:),nTrlPerCmp(i));
        end

        PC=PsyCurveData.percentCorrect(stdX,cmpX,RDat);
        out=PsyCurveData(stdX,cmpX,RDat);
        out.G=g;
    end
    function [S,G]=gen(stdX,mu,sigma,cr,nTrlPerCmp)
        [S,G]=PsyCurveData.genNorm(mu,sigma,cr,nTrlPerCmp);

        sz=size(S.RCmpChs);

        nS=numel(stdX);

        cunq=unique(S.cC);
        nC=numel(cunq);

        S.stdX=zeros(sz);
        S.cmpX=zeros(sz);
        for is = 1:nS
            indS=S.sC==is;
            S.stdX(indS)=stdX(is);
            for ic = 1:nC
                indC=S.cC==ic;
                S.cmpX(indS & indC)=G.mu(is,ic)+stdX(is);
            end
        end

    end
    function [S,G]=genNorm(mu,sigma,cr,nTrlPerCmp)
        % everything has been normalized
        % cr = 0 = unbiased
        % mu - relative to center
        %
        % mu - at least a row vec
        % sigma - at least a col vec

        % Correct dimensions
        if Vec.isCol(mu)
            mu=mu';
        end
        if Vec.isRow(sigma)
            sigma=sigma';
        end

        nCmp=size(mu,2);
        nStd=size(sigma,1);
        n=nCmp*nStd;

        % correct sizes
        if size(mu,1)==1
            mu=repmat(mu,nStd,1);
        end
        if size(sigma,2)==1
            sigma=repmat(sigma,1,nCmp);
        end
        if numel(cr)==1
            cr=repmat(cr,nStd,nCmp);
        elseif size(cr,2)==1
            cr=repmat(cr,1,nCmp);
        elseif size(cr,1)==1
            cr=repmat(cr,nStd,1);
        end
        if numel(nTrlPerCmp)==1
            nTrlPerCmp=repmat(nTrlPerCmp,nStd,nCmp);
        elseif size(nTrlPerCmp,2)==1
            nTrlPerCmp=repmat(nTrlPerCmp,1,nCmp);
        elseif size(nTrlPerCmp,1)==1
            nTrlPerCmp=repmat(nTrlPerCmp,nStd,1);
        end

        DVDat=[];
        RDat=[];
        inds=[];
        for i = 1:n
            % DV DATA
            dvdat=normrnd(mu(i),sigma(i),nTrlPerCmp(i),1);
            DVDat=[DVDat; dvdat];

            % R DATA
            rdat =bsxfun(@gt,dvdat,cr(i));
            RDat=[RDat; rdat];

            [s1,s2]=ind2sub([nStd nCmp],i);
            % stdX cmpX
            inds=[inds; repmat([s1 s2],nTrlPerCmp(i),1)];
        end
        S.RCmpChs=RDat;
        S.sC=inds(:,1);
        S.cC=inds(:,2);

        G.mu=mu;
        G.sigma=sigma;
        G.cr=cr;
        G.DVDat=DVDat;

    end
    function [PC,N1,N,N0]=percentCorrect(stdX,cmpX,RCmpChs)
        % STANDARD VALUES
        XstdUnq = unique(stdX);

        % LOOP OVER STANDARDS
        for s = 1:length(XstdUnq)
            % INDICES FOR EACH STANDARD
            indS = stdX == XstdUnq(s);
            % COMPARISON VALUES
            XcmpUnq(:,s) = unique(cmpX(indS));
            % LOOP OVER COMPARISONS
            for c = 1:length(XcmpUnq(:,s))
                % INDICES IN STD / CMP CONDITION
                indCnd  = stdX ==XstdUnq(s) & cmpX==XcmpUnq(c,s);
                % TOTAL NUMBER OF TRIALS IN CONDITION
                N(c,s)  = sum( indCnd );
                % TOTAL NUMBER OF CMP CHOSEN IN CONDITION
                N1(c,s) = sum( RCmpChs(indCnd)==1 );
                % TOTAL NUMBER OF STD CHOSEN IN CONDITION
                N0(c,s) = sum( RCmpChs(indCnd)==0 );
            end
        end

        PC = N1./N;
    end
end
end
