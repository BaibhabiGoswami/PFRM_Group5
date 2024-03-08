function [BSretMeans,BSvcov] = BS (RetData,cov,retMeans)
% Using Bayes-Stein estimator

% Weights Hist MVP
[n,col]=size(RetData)
cov_inv = inv(cov);
ones_vector = ones(col,1)
wh = cov_inv * ones_vector/(ones_vector'* cov_inv * ones_vector)
% Volatility of Hist MVP
vh = (wh' * cov) * wh
% Mean excess-return of Hist MVP
rh = wh' * retMeans'
% T (number of observations)= n
% N (number of assets)= col
% Lambda
lambda =((col+2)*(n+1))/(((retMeans-(rh*ones_vector'))*cov_inv)*(retMeans-(rh*ones_vector'))'*(n-col-2))
% Psi (= Lambda/(T+Lambda)))
Psi = lambda / (n+lambda)
% Expected excess-returns Bayes-Stein
BSretMeans = (1-Psi) * retMeans + Psi * rh * ones'
% Covariance Matrix Bayes-Stein
variable1 =lambda/(n*(n+1+lambda))
variable2 = (ones_vector*ones_vector')/(ones_vector'*cov_inv*ones_vector)
BSvcov=cov*(1+(1/(n+lambda)))+variable1*variable2