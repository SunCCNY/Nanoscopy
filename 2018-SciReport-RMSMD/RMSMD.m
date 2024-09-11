% [D D2]=RMSMD(S,X) 
%
% Calculate RMSMD - Root mean square minimum distance of two point sets S and X
% 
% For technical details, see:
% [1] Sun, Y. Mean Square Minimum Distance: a Quality Metric for Localization 
% Nanoscopy Images. Submitted to Scientific Reports (2018). 
% 
% Input:
%   S         - size(S)=[dim,NS], NS dim-dimensional vectors, each represents a point
%   X         - size(X)=[dim,NX], NX dim-dimensional vectors, 
%
% Output:
%   D         - RMSMD between S and X 
%   D2        - MSMD between S and X, that is, D^2
%
% 04/08/2017
% 08/19/2019 Same as RMSMD_Core.m 
% Yi Sun

function [D,D2]=RMSMD(S,X) 

[dim,M]=size(S) ; % dim - dimension of S, NS - # of points in S
[din,K]=size(X) ; % din - dimension of X, NX - # of points in X
if dim~=din
  fprintf(1,'S and X have different dimensions.\n')
  return ;
end
if M<=0||K<=0
  fprintf(1,'S and X cannot be empty.\n') ;
  return ;
end
%% Search minimum-distance points 
% For each x in X, search s in S with minimum distance 
% S_=zeros(size(X)) ; % points in S with minimum distance for each xi 
D_X2S=zeros(1,K) ; % minimum distance for each xi to S 
for i=1:K
 tmp=S-X(:,i)*ones(1,M) ; 
 d=sum(tmp.^2) ;
 [D_X2S(i),j]=min(d) ; 
% S_(:,i)=S(:,j) ; 
end
% For each s in S, search x in X with minimum distance 
% X_=zeros(size(S)) ; % points in X with minimum distance for each si 
D_S2X=zeros(1,M) ; % minimum distance for each si to X 
for i=1:M
 tmp=X-S(:,i)*ones(1,K) ;
 d=sum(tmp.^2) ;
 [D_S2X(i),j]=min(d) ; 
% X_(:,i)=X(:,j) ; 
end
D2=(sum(D_S2X)+sum(D_X2S))/(M+K) ;
D=sqrt(D2) ; 

end
