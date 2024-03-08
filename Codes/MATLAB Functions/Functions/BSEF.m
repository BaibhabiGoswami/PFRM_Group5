function [meanGMVP,stddevGMVP,meanTP,stdDevTP,varpor_r,retport_r,retMeansTarget,filteredRet,filteredVar] = BSEF (xretMeans,xmvarcov,up,low)
%riskfree is timeseries
% xret
%RF=RiskFree
%FactorsFive=FactorsFivedata(3:722,2:6)
%xret=RiskyReturns-RF
[T,N] = size(xretMeans) ; 

% 5 Factors expected returns and excess returns
%retMeans  = mean(RiskyReturns) ;
%xretMeans= mean( xret ) ; % I will need it when calculating weights for the tangency portfolio

% Standard deviation of returns
%StdDev = std( retMeans ) ;

% demeaned returns and excess returns for calculating variances and covariances
%demeaned_ret= RiskyReturns - retMeans; 
%demeaned_xret = xret - xretMeans; 

% Variance-covariance matrix of returns and excess returns
%mvarcov = cov( demeaned_ret ); 
%xmvarcov = cov( demeaned_xret ) ; % I will need it when calculating weights for the tangency portfolio

toc
%% GMVP

% weights of GMVP calculated using the analytical formula: WGmvp = ((V^(-1))1)/(1'(V^(-1))1)
%weights unscaled
weightsGMVP = (xmvarcov^(-1))*ones(N,1)  ;

%weights scaled to sum to 1
weightsGMVP = weightsGMVP./sum(weightsGMVP) ;

% GMVP portfolio statistics - expected retun, variance and volatility
meanGMVP = xretMeans*weightsGMVP ;  
varGMVP  = weightsGMVP'*xmvarcov*weightsGMVP ;
stddevGMVP  = sqrt(varGMVP) ;

toc

% analytical solution for weights of tangency portfolio 

%weights unscaled
weightsTP = xmvarcov^(-1)*xretMeans';

%weights scaled to sum to 1
weightsTP = weightsTP./sum(weightsTP);

% Tangency portfolio statistics
meanTP = xretMeans*weightsTP ;
varTP  = weightsTP'*xmvarcov*weightsTP ;
stdDevTP  = sqrt(varTP) ;

%% Risky EF
%low=min(xret)
step = (up-low)/100
retMeansTarget = (low:step:up)

mm=length(retMeansTarget)
    varpor=zeros(1,mm);
    weights = zeros(N,mm);
    c = zeros(N,1);
    Aeq_r=[xretMeans;ones(1,N)];
    filteredRet = [];
    filteredStdDev = [];
    filteredVar=[];
    %risky assets MVF
    for i=1:mm
        beq_r = [retMeansTarget(i);1];
        vlb = -ones(N,1)*10000
        vub = ones(N,1)*10000;
        x0 = repmat(1/N, N, 1);
        options = []; 
        options = optimset('display', 'on','largescale', 'off');
        x_r = quadprog(xmvarcov,c,[],[],Aeq_r,beq_r,vlb,vub,x0,options);
        varpor_r(i) = x_r'*xmvarcov*x_r;
        weights_r(:,i) = x_r;
    end
    retport_r=xretMeans*weights_r
    stdDevpor_r=sqrt(varpor_r)

    [minStdDev, minIndex] = min(stdDevpor_r);
    minStdDevReturn = retport_r(minIndex);

    % Filter for portfolios with expected returns >= minStdDevReturn
    for i = 1:mm
        if retport_r(i) >= minStdDevReturn
            filteredRet = [filteredRet, retport_r(i)];
            filteredStdDev = [filteredStdDev, stdDevpor_r(i)];
        end
    end
    filteredVar=filteredStdDev.^2