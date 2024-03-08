function [rejected, pval, testStat] = robustSharpe(data,alpha,H0,M,bl,kernel,extsim,bootMat)
% [rejected, pval, teststat] = robustSharpe(data,alpha,H0,M,bl,kernel,extsim,bootMat)
%  minimal input             = robustSharpe(data)
% Tests whether the difference of two Sharpe ratios is statistically significantly 
% different from H0. Tested as two-sided hypothesis by simulating M datasets by the 
% circular block bootstrap, taking critical values as the empirical quantiles of the 
% simulated datasets.
% The Quadratic Spectral Kernel is used in the prewhitened kernel variance 
% estimator. Based on Ledoit&Wolf (2008).
%
% Inputs:
%   data:       [Tx2] matrix of excess returns in format 0.02 equiv. 2% excess return
%               (if you feed in 1.02 equiv. 2% excess return, the internally computed
%               Sharpe ratios will be wrong)
%   alpha:      fixed significance level; default value = 0.05
%   H0:         null hypothesized value for the value of Sharpe ratios difference;
%               default value = 0 
%   M:          number of bootstrap iterations; default value = 5,000
%   bl:         block size in Circular Block Bootstrap.
%               Use routine optimalblrobustSharpe.m to determine optimal block size.
%				If no block size is specified, optimalblrobustSharpe.m is called automatically,
%				with default candidate block sizes 1,3,6,10,15
%   kernel:     'G' for Gallant/Parzen, 'QS' for Quadratic spectral (default)
%				Note: The R code of the paper uses the Gallant/Parzen kernel.
%   extsim:     1 if the indices matrix bootMat in the circluar block bootstrap is fed
%               in rather than simulated in robustSharpe itself, 0 else
%               useful to achieve comparability of results based on other implementations
%   bootMat:    exogenous indices matrix in circular block bootstrap of size [MxT] or 0
%               where M is number of CBB iterations, T is time series length
%
% Outputs:
%   rejected:   1 if H0 was rejected at significance level alpha, 0 else
%   pval:       p-value 
%   teststat    test statistic
%
% ©2009 Dan Wunderli, Institute for Empirical Research in Economics, U Zurich

%Defaults
if not(ismember('stud',who)), stud=1; end
if not(ismember('M',who)), M=5000; end
if not(ismember('bl',who)), fprintf('%s \n','Computation of optimal block size:'), bl=optimalblrobustSharpe(data,0,2000,200); fprintf('%s \n','Testing of H0:'), end
if not(ismember('alpha',who)), alpha=0.05; end
if not(ismember('H0',who)), H0=0; end
if not(ismember('extsim',who)), extsim=0; end
if not(ismember('bootMat',who)), bootMat=0; end
if not(ismember('kernel',who)), kernel='QS'; end

%Start
S=size(data,2); t=size(data,1);

w_T=0;
%Computation of studentized test statistic and Generation Circular Block Bootstrap Index Matrix
X=data;

theta_0 = ones(1,S)*H0;
mu = mean(X,1);
mui = mean(X(:,1),1)/std(X(:,1))-mean(X(:,2),1)/std(X(:,2));

%%%needed input arguments for andmon implementation of Mike Cliff, V1.1
gmmopt.prt             = 0;
gmmopt.aminfo.p        = 1;
gmmopt.aminfo.q        = 0;
gmmopt.aminfo.vardum   = 0;
gmmopt.aminfo.kernel   = kernel;
gmmopt.aminfo.nowhite  = 0;
gmmopt.aminfo.diagdum  = 0;
gmmopt.plot            = 0;

nu = [mean(X(:,1)); mean(X(:,2)); mean(X(:,1).^2); mean(X(:,2).^2)];
nux = [X(:,1)-nu(1) X(:,2)-nu(2) X(:,1).^2-nu(3) X(:,2).^2-nu(4)];
nabla = [nu(3)/((nu(3)-nu(1)^2)^1.5) -nu(4)/((nu(4)-nu(2)^2)^1.5) -0.5*nu(1)/((nu(3)-nu(1)^2)^1.5) 0.5*nu(2)/((nu(4)-nu(2)^2)^1.5)]';
%%%Prewhiten data with VAR(1) model, estimate HAC kernel estimator
%%%using AR(1) models as univariate approximating parametric models
dataVAR = vare2(nux,1,0);          %alternatively by vare(nux,1)
THETA = horzcat(dataVAR(1:4).beta)';
[U,S,V]=svd(THETA);
for i=1:size(THETA,2)
    if S(i,i)>0.97
        S(i,i)=0.97;
    elseif S(i,i)<-0.97
        S(i,i)=-0.97;
    end
end
THETA = U*S*V';

u = (nux(1+1:end,:)' - THETA*nux(1:end-1,:)')'; %equals V.star in R implementation

[band,Z,alphaa,num,den,bandoverall,thet,sigma2] = andmon6cvm2(gmmopt,u,THETA);    %covariance matrix by kernel-based HAC estimator of Andrews-Monahan (1992). andmon6 is an altered version of Mike Cliff's andmon function Version 1.1
varde = (nabla'*((eye(size(u,2))-THETA)\Z/(eye(size(u,2))-THETA)')*nabla);

z_T = abs(mui)./(sqrt(varde/t));        %studentization of `raw' test statistic
mu;										%means of two return time series
mui;									%Difference of Sharpe ratios
sigi = sqrt(varde)./sqrt(t);            %HAC std estimate
r = mod(t,bl); L = floor(t/bl);
Xbootind = zeros(M,t);

if extsim==0
    for m = 1:M                         %generation of M cbb matrices X_T*m, 1<=m<=M
        [ind,L,r] = cbb_seq(t,bl);
        Xbootind(m,:) = ind;
    end
elseif extsim==1
    Xbootind=bootMat;
else error('extsim=1 or extsim=0 accepted')
end

%bsstat is a matrix with corresponding studentized test statistics for each bootstrap iteration (row)
%muboot contains CBB simulated excess returns of two assets,
%sigboot is HAC std estimate of difference of two Sharpe ratios
[bsstat,muboot,sigboot] = bsstats11(X,Xbootind,mui,bl,L,r,t);

[d_K,rejected] = critvalues1(bsstat,mui,sigi,theta_0,alpha);
%computes critical value and tests H0

testStat=z_T;
pval=(length(find(bsstat>=abs(z_T)))+1)/(M+1);
pValue = pval;

fprintf('%s \n',horzcat('H0: ',num2str(H0),', stud: ',num2str(stud),', M: ',num2str(M),', bl: ',num2str(bl),', alpha: ',num2str(alpha)))

fprintf('%s \n',horzcat('test statistic: ',num2str(testStat,4), ...
		    ', p-value: ',num2str(pval,3),', rejected: ',num2str(rejected)))