function [d_K,rej] = critvalues1(bsstat,mu,sigi,theta_0,alph) %bsstat,S,mu,sigi,theta_0,k,alph,Nmaxh
% iteratively calculates the critical value d_K from the bootstrap matrix of studentized test
% statistics called bsstat (rows are simulation runs, colums variables). Does iterative testing. S is size(data,2), sigi are estimated
% s.d. in original word, theta_0 are hypothesized values, k, alph, Nmax as in kStepM.

zma = bsstat;

zmat = zma;              
d_K = quantileR(zmat,1-alph,1);   

%invert generalized confidence regions, i.e. reject hypothesis s, 1<=s<=S
rej = zeros(1);
    if  d_K < abs(mu-theta_0)/sigi
        rej = 1;
    end

rej;