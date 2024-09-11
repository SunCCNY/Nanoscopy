% [Xp,Ki,Ip]=partitionX(S,X,Kai) 
%
% Partition a set of estimated locations according to a set of truth
% locations S and number of activations in a data movie
% 
% For technical details, see:
% [1] Sun, Y. Mean Square Minimum Distance: a Quality Metric for Localization 
% Nanoscopy Images. Submitted to Scientific Reports (2018). 
% 
% Input:
%   S     - size(S)=[D,M]; M D-dimensional true emitter locations 
%   X     - size(X)=[D,K]; K D-dimensional estimated locations
%   Kai   - size(1,M); Kai(i) is number of activations for ith emitter 
%           in a data movie
%
% Output:
%   Xp    - size(D,N) with K=sum(Ki); Xp(:,1:K1) are estimated locations 
%           that are assigned to 1st emitter after partition, and rest are 
%           stacked in the order of emitter indices
%   Ki    - size(1,M), Ki(m) is number of estimated locations for partition 
%           with mth emitter
%   Ip    - size(2,Kp) with Kp=sum(Ki), ith location in Xp is Ip(1,i)th 
%           location in X and is assinged to Ip(2,i)th emitter. Ip is used 
%           for computation of Precision, Recall, F1
% 
% Yi Sun
% 09/10/2021, 02/14/2022

function [Xp,Ki,Ip]=partitionX(S,X,Kai) 

[Ds,M]=size(S) ;    % ds = dimension of S, M = # of locations in S
[Dx,K]=size(X) ;    % dx = dimension of X, K = # of locations in X
if Ds~=Dx
  fprintf(1,'S and X have different dimensions.\n')
  return ;
end
if M<=0||K<=0
  fprintf(1,'S and X cannot be empty.\n') ;
  return ;
end
D=zeros(M,K) ;      % matrix of square distances between S(:,i) and X(:,j) 
for m=1:M
  tmp=X-S(:,m)*ones(1,K) ; 
  D(m,:)=sum(tmp.^2,1) ; 
end
Db=D ;              % backup of D
if K<2*M            % a joint moive localization (JML) algorithm
  beta=ones(1,M)/M ;
else                % a frame by frame localizaation (FFL) algorithm
  beta=Kai/sum(Kai) ;
end
Ki=round(K*beta) ;
Xt=zeros(Ds,max(Ki),M) ; 
pt=zeros(size(Kai)) ;  % pt(m) is a pointer to number of locations in Xt(:,:,m) 
It=zeros(2,K+M) ;   % ith location in Xp is It(i)th location in X
TKi=sum(Ki) ;       % total # of locations supposed to partition
if K<=TKi
  for n=1:K
    [V,r]=min(D) ;
    [~,c]=min(V) ;
    m=r(c) ;        % cth location has minimum distance to mth emitter
    pt(m)=pt(m)+1 ;
    Xt(:,pt(m),m)=X(:,c) ;
    D(:,c)=Inf ;    % remove cth location
    if pt(m)==Ki(m) % number of locations for mth emitter reaches maximum
      D(m,:)=Inf ;  % remove mth emitter
    end
    It(1,n)=c ;
    It(2,n)=m ;
  end
end
if K>TKi
  Dt=D ; 
  for n=1:TKi
    [V,r]=min(D) ;
    [~,c]=min(V) ;
    m=r(c) ;        % cth location has minimum distance to mth emitter
    pt(m)=pt(m)+1 ;
    Xt(:,pt(m),m)=X(:,c) ;
    D(:,c)=Inf ;    % remove cth location
    if pt(m)==Ki(m) % number of locations for mth emitter reaches maximum
      D(m,:)=Inf ;  % remove mth emitter
    end
    It(1,n)=c ;
    It(2,n)=m ;
    Dt(:,c)=Inf ;   % remove only cth location but not emitter from Dt 
  end
  % Note: now D(i,j)=Inf for all i,j, and Dt(i,j)=Inf only for partitioned locations
  D=Dt ;            % partition remainder
  for n=TKi+1:K
    [V,r]=min(D) ;
    [~,c]=min(V) ;
    m=r(c) ;        % cth location has minimum distance to mth emitter
    pt(m)=pt(m)+1 ;
    Xt(:,pt(m),m)=X(:,c) ;
    D(:,c)=Inf ;    % remove cth estimated location
    D(m,:)=Inf ;    % remove mth emitter anyway. Note: K-Tki<=M
    It(1,n)=c ;
    It(2,n)=m ;
  end
end
% If no estimate is assigned to mth emitter, the closest estimate is 
% assigned to mth emitter -> the estimate is equivalent to two identical
% estimates
for m=1:M
  if pt(m)==0 
    [~,c]=min(Db(m,:)) ;
    pt(m)=pt(m)+1 ;
    Xt(:,pt(m),m)=X(:,c) ;
    n=n+1 ;         % n continues with above loop
    It(1,n)=c ; 
    It(2,n)=m ; 
  end
end
Ki=pt ; 
Kp=sum(Ki) ;        % Kp might be greater than K due to the last step above
Xp=zeros(Ds,Kp) ; 
p=0 ;
for m=1:M
  Xp(:,p+1:p+Ki(m))=Xt(:,1:Ki(m),m) ;
  p=p+Ki(m) ;
end
Ip=It(:,1:n) ;      % Note: n=Kp

end
