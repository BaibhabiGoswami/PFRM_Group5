function [optbl blcand errr] = optimalblrobustSharpe(data,silent,M,MBM,blcand,avbl,alpha,H0,Tstart,extsimSB,extSB,extsimMB,extMB)
% [optbl blcand errr] = optimalblrobustSharpe(data,silent,M,MBM,blcand,avbl,alpha,H0,Tstart,extsimSB,extSB,extsimMB,extMB)
% minimal input       = optimalblrobustSharpe(data)
% In robustSharpe, confidence regions are simulated. To find the optimal block size, 
% one first simulates datasets by a parametric VAR(1) bootstrap with stationary bootstrapping
% the residuals. Then one constructs confidence regions by circular block bootstrap
% constructed critical values and observes how often rejections of the null hypothesis occur.
% The optimal block size is the one that minimizes the deviation of the simulated
% coverage from the desired coverage 1-alphaa i.e. minimizes the deviation of the
% desired error rate alphaa from the simulated error rate errr.
%
% Inputs:
%   data:   [Tx2] matrix of excess returns
%   silent: 1 for yes, 0 for no. Default is not silent, so that
%           estimated time left in seconds and current error rates are printed
%   M:      number of `outer' bootstrap iterations (stationary bootstrap)
%           default = 2,000;
%   MBM:    number of `inner' bootstrap iterations (circular blocks bootstrap);
%           default = 200
%   blcand: candidate block sizes, default = [1 3 6 10 15]
%   avbl:   average blocksize in first order bootstrap (stationary bootstrap)
%   alpha:   fixed significance level; default = 0.05
%   H0:     null hypothesized value, typically zero in this context
%   Tstart: in the parametric bootstrap, e.g. Tstart=1000 to make starting values irrelevant
%   extsimSB:   1 if the indices matrix extSB in the stationary bootstrap is fed
%               in rather than simulated in optimalblrobustSharpe itself, 0 else
%   extSB:      exogenous indices matrix in stationary bootstrap of size [M x T] or 0
%   extsimMB:   1 if the indices matrix extSB in the circular blocks bootstrap is feeded
%               in rather than simulated in optimalrobustSharpe itself, 0 else
%   extMB:      exogenous indices matrix in circular block bootstrap of size [MBM x T] or 0
%
% Outputs:
% (during execution, estimated time left in seconds and current error rates are printed)
%   optbl:  the optimal blocksize so that the simulated coverage is closest to the desired one
%   blcand: blcand vector that was put into the algorithm
%   errr:   error rate achieved in the simulation, for each candidate blocksize
%
% ©2009 Dan Wunderli, Institute for Empirical Research in Economics, U Zurich

format short;

if not(ismember('MBM',who)), MBM=200; end
if not(ismember('M',who)), M=2000; end
if not(ismember('alpha',who)), alpha=0.05; end
if not(ismember('avbl',who)), avbl=5; end
if not(ismember('blcand',who)), blcand=[1 3 6 10 15]; end
if not(ismember('H0',who)), H0=0; end
if not(ismember('Tstart',who)), Tstart=1000; end
if not(ismember('extsimSB',who)), extsimSB=0; end
if not(ismember('extSB',who)), extSB=0; end
if not(ismember('extsimMB',who)), extsimMB=0; end
if not(ismember('extMB',who)), extMB=0; end

varlag=1;
t = size(data,1)-varlag; T = size(data,1); S = size(data,2);
% datac = data - repmat(mean(data,1),[T 1]);

theta_0star = mean(data(:,1))/std(data(:,1)) - mean(data(:,2))/std(data(:,2));

Pt = vare2(data,1,1);                        %!!!!varsvdfit numerically more advisable here!!!

yhat=horzcat(Pt(1:S).yhat);
A=horzcat(Pt(1:S).beta)';
% A = A(:,1:end-1);                               %cuts off estimate of constant

resid = data(1+varlag:end,:) - yhat;
residc = resid;
% residc = resid - repmat(mean(resid,1),[t 1]);   %!!!better numerical centering?

% indmatSB = zeros(M,t+Tstart);

%% needed input arguments for andmon implementation of Mike Cliff, V1.1
gmmopt.prt             = 0;
gmmopt.aminfo.p        = 1;
gmmopt.aminfo.q        = 0;
gmmopt.aminfo.vardum   = 0;
gmmopt.aminfo.kernel   = 'G';
gmmopt.aminfo.nowhite  = 0;
gmmopt.aminfo.diagdum  = 0;
gmmopt.plot            = 0;

%%%%%%%%%%% Simulation of data sets by Stationary Bootstrap
Out = zeros(1,length(blcand));
tic
for m=1:M       %%% first level bootstrap
    timeb=toc;
    if extsimSB==0
        sbb_seq = sbbseq(t,Tstart,avbl);        % stat. BB sequence of length T+Tstart
        % indmatSB(m,:) = sbb_seq;
    else sbb_seq=extSB(m,:);
    end
    residstar = residc(sbb_seq,:);

    xold = zeros(S,1); xboot = zeros(t+Tstart,S);

    for j=1:(t+Tstart)
        xnew = A(:,3)+ A(:,[1 2])*xold + residstar(j,:)';
        xboot(j,:) = xnew';
        xold = xnew;
    end
    xboot = xboot(Tstart:t+Tstart,:);           % throw away first Tstart-1 values

    %%% Studentization of simulated data if stud==1


    % theta_0 = ones(1,S)*H0;       % vector of hypothesized values w.r.t. test statistic sorted IX!!!
    X = xboot;
    mub = mean(X,1);
    mubi = mean(X(:,1))/std(X(:,1)) - mean(X(:,2))/std(X(:,2));
        vard = zeros(1,S);

       
            nu = [mean(X(:,1)); mean(X(:,2)); mean(X(:,1).^2); mean(X(:,2).^2)];
            nux = [X(:,1)-nu(1) X(:,2)-nu(2) X(:,1).^2-nu(3) X(:,2).^2-nu(4)];
            nabla = [nu(3)/((nu(3)-nu(1)^2)^1.5) -nu(4)/((nu(4)-nu(2)^2)^1.5) -0.5*nu(1)/((nu(3)-nu(1)^2)^1.5) 0.5*nu(2)/((nu(4)-nu(2)^2)^1.5)]';
            %%% Prewhiten data with VAR(1) model, estimate HAC kernel estimator
            %%% using AR(1) models as univariate approximating parametric models
            dataVAR = vare2(nux,1,0);          % alternatively by vare(nux,1)
            THETA = horzcat(dataVAR(1:4).beta)';
            [U,Q,V]=svd(THETA);
            for i=1:size(THETA,2)
                if Q(i,i)>0.97
                    Q(i,i)=0.97;
                elseif Q(i,i)<-0.97
                    Q(i,i)=-0.97;
                end
            end
            THETA = U*Q*V';
            

            u = (nux(1+1:end,:)' - THETA*nux(1:end-1,:)')'; % equals V.star in R implementation

            [band,Z,alphaa,num,den,bandoverall,thet,sigma2] = andmon6cvm2(gmmopt,u,THETA);    % covariance matrix by kernel-based HAC estimator of Andrews-Monahan (1992). andmon6 is an altered version of Mike Cliff's andmon function
            % Version 1.1
            % recoloring
            vard = (nabla'*((eye(size(u,2))-THETA)\Z/(eye(size(u,2))-THETA)')*nabla);

        sigbl = sqrt(vard./T);
    
    
    for b=1:length(blcand)
        blc=blcand(b);

        % second-level bootstrap
        Xbootind=zeros(MBM,T);
        for i=1:MBM
            if extsimMB==0
                [ind,L,r] = cbb_seq(T,blc);
                Xbootind(i,:) = ind;
            elseif extsimMB==1
                Xbootind(i,:) = extMB(i,:);
            end
        end

        % computation of simulated coverage and ghat(b)
            [bsstatistics,muboot,sigboot] = bsstats11(X,Xbootind,mubi,blc,floor(T/blc),mod(T,blc),T);

            [d_1,rej] = critvalues1(bsstatistics,mubi,sigbl,theta_0star,alpha);

            Out(b) = Out(b) + rej;
    end
    timee=toc;
    NRP = sum(Out,1)./m;
    time=['Estimated time remaining: ' num2str((timee-timeb)*(M-(m-1)),3) ...
	  's, null rejection probs: ' num2str(NRP,3)];
    fprintf('%s \n',time)
end

NRP = sum(Out,1)./M

optbl = blcand(find(abs(NRP-alpha)==min(abs(NRP-alpha))))

csvwrite('optimalblrobustSharpe.csv',[blcand NRP optbl]);