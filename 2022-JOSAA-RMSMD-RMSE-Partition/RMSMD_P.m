% D=RMSMD_P(S,X,Ki) 
%
% Calculate RMSMD for two point sets S and X with known partition of X
% 
% For technical details, see:
% [1] Sun, Y. Mean Square Minimum Distance: a Quality Metric for Localization 
% Nanoscopy Images. Submitted to Scientific Reports (2018). 
% 
% Input:
%   S     - size(S)=[D,M]; M D-dimensional true emitter locations 
%   X     - size(X)=[D,K]; K D-dimensional lcoation estimates. X(:,1:Ki(1))
%           are estimates for 1st emitter. Other estimates are stocked
%           accordingly.
%   Ki    - size(1,M); Ki(m) is number of estimates for mth emitter in X
%
% Output:
%   D         - RMSMD_P between S and X 
%
% 09/14/2021
% Yi Sun

function D=RMSMD_P(S,X,Ki) 

[ds,M]=size(S) ; % ds - dimension of S, M - # of points in S
[dx,K]=size(X) ; % dx - dimension of X, K - # of points in X
if ds~=dx
  fprintf(1,'S and X have different dimensions.\n')
  return ;
end
if M<=0||K<=0
  fprintf(1,'S and X cannot be empty.\n') ;
  return ;
end
%% Search minimum-distance points 
% For each x in X, search s in S with minimum distance according to
% partition Xi
D_X2S=0 ; % minimum distance for each xi to S 
p=0 ; 
for m=1:M
  D_X2S=D_X2S+sum(sum((X(:,p+1:p+Ki(m))-S(:,m)*ones(1,Ki(m))).^2)) ;
  p=p+Ki(m) ;
end
% For each s in S, search x in Xi with minimum distance 
D_S2X=0 ; % minimum distance for each si to X 
p=0 ; 
for m=1:M
  d=sum((X(:,p+1:p+Ki(m))-S(:,m)*ones(1,Ki(m))).^2,1) ;
  D_S2X=D_S2X+min(d) ;
  p=p+Ki(m) ;
end
D2=(D_S2X+D_X2S)/(K+M) ;
D=sqrt(D2) ; 

end
