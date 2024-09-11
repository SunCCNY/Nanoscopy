% h=RMSE_P(S,X,Ki) 
%
% Calculation of RMSE with known partition of estimates.
% 
% For technical details, see:
% [1] Y. Sun. Potential quality improvement of stochastic optical 
% localization nanoscopy images obtained by frame by frame localization 
% algorithms. Scientific Reports, vol. 10, no. 1, pp. 11824, July 16, 2020.
% 
% 
% Input:
%   S     - size(S)=[D,M]; M D-dimensional true emitter locations 
%   X     - size(X)=[D,K]; K D-dimensional lcoation estimates. X(:,1:Ki(1))
%           are estimates for 1st emitter. Other estimates are stocked
%           in order
%   Ki    - size(1,M); Ki(m) is number of estimates for mth emitter in X
%
% Output:
%   h    - RMSE with known partition of X
%
% Yi Sun
% 09/14/2021

function h=RMSE_P(S,X,Ki)

[Ds,M]=size(S) ; % ds = dimension of S, M = # of locations in S
[Dx,K]=size(X) ; % dx = dimension of X, K = # of locations in X
if Ds~=Dx
  fprintf(1,'S and X have different dimensions.\n')
  return ;
end
if M<=0||K<=0
  fprintf(1,'S and X cannot be empty.\n') ;
  return ;
end
h=0 ; 
p=0 ;
for m=1:M
  h=h+sum(sum((X(:,p+1:p+Ki(m))-S(:,m)*ones(1,Ki(m))).^2)) ;
  p=p+Ki(m) ;
end
h=sqrt(h/K) ;

end
