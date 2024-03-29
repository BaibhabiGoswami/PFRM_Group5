function  q = quantileR(x,p,method)
% PURPOSE: compute empirical quantile (percentile). Default
% here is equivalent to default of function `quantile' in R (project-r.org).
% ----------------------------------------------------
% USAGE:   q = quantile(x,p,method)
% where:   x = matrix or vector 
%          p = percent
%     method = 1,2,3
%    1. Interpolation so that F(X_(k)) == (k-1)/(n-1) (default)
%    2. Interpolation so that F(X_(k)) == k/(n+1).
%    3. The least number q such that at least a part p of x 
%      is less than or equal to q. 
% -----------------------------------------------------
% RETURNS: q = empirical quantile
% -----------------------------------------------------
% NOTES:  If input x is a matrix then the quantile is computed for 
%        every column. Input p may be vector also.     

if nargin<3, method=1; end
if min(size(x)) == 1
   x = x(:);
   q = zeros(size(p));
else
   q = zeros(length(p),size(x,2));
end
if min(size(p)) > 1
   error('Not matrix p input')
end
if any(p>1|p<0)
   error('Input p is not between 0 and 1')
end

x = sort(x); 
p = p(:);
n = size(x,1);
if method == 3
   qq = x(ceil(min(max(1,p*n),n)),:);
else                         % Method 1 is from Hjort's "Computer
   x = [x(1,:); x; x(n,:)];  % intensive statistical methods" page 102
   if method == 2
      i = p*(n+1)+1;
   else
      i = p*(n-1)+2;
   end
   iu = ceil(i);
   il = floor(i);
   d = (i-il)*ones(1,size(x,2));
   qq = x(il,:).*(1-d)+x(iu,:).*d;
end

q(:) = qq;

