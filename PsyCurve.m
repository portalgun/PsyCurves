classdef PsyCurve < handle
properties
end
methods(Static)
    function [PC,DP] = genGauss(X,mu,sigma,bet,nIntrvl)
        DP = sign(X-mu).*(abs(X-mu)./sigma).^bet; % abs() prevent complex numbers w. some betas
        PC = normcdf( 0.5.*sqrt(nIntrvl).*DP,0,1); % sign() reinstates sign of Xcmp-muFit
    end
    function T=sigma2thresh(sigma,dpcrit,bet)
        if nargin < 2
            dpcrit=1;
        end
        if nargin < 3
            bet=1;
        end

        % XXX T  = sigma.*DPcrt.^(1./bet); from psygengauss
        T=sqrt(sigma.^2.*dpcrit.^(1./bet));
        % WHAT I HAD
        %T=sqrt(sigma.^2./(dpcrit.^(1./bet)));
    end
    function sigma=thresh2sigma(T,dpcrit,bet)
        if nargin < 2
            dpcrit=1;
        end
        if nargin < 3
            bet=1;
        end

        sigma=sqrt(T.^2./(dpcrit.^(1./bet)));
        % WHAT I HAD
        %sigma=sqrt(thresh.^2.*dpcrit.^(1./bet));

    end
    function [sigma,t]=lin2sigma(x,m,y,bet,bLog,dpcrit);
        if nargin < 4
            bet=1;
        end
        if nargin < 5
            bLog=true;
        end
        if nargin < 6
            dpcrit=1;
        end

        %Y=log(obj.T(:,m,b,k,r));
        %fun=@(mb) mean(abs(ffun(mb) - y).^2);
        %TFit(:,m,b,k,r)=exp(ffun(fout));
        if bLog
            t=exp(m.*x + y);
            %t=exp(m.*x) + y;
            %t=m.*x + y;
        else
            t=m .* x + y;
        end
        sigma=PsyCurve.thresh2sigma(t,dpcrit,bet);
        %sigma=t;
    end
end
methods
    function obj=PsyCurve()
    end
end
end
