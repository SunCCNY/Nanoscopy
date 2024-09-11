% [F1,Pre,Rec]=MCmetrics(Ip,Ix)
%
%   Compute metrics for multiclass classification. Assume that the original 
% samples in X is classified and arranged in Xp.
% 
% Input:
%   Ix    - size(1,K), ith sample in X belongs to Ix(i)th class
%   Ip    - size(2,Kp), ith sample in Xp is Ip(1,i)th sample of X and is 
%           classified into Ip(2,i)th class
%
% Output:
%   F1    - F1 score over all classes
%   Pre   - Precision averaged over all classes
%   Rec   - Recall averaged over all classes
%
% Yi Sun
% 02/14/2022

function [F1,Pre,Rec]=MCmetrics(Ip,Ix)

M=max(Ix) ;       % # of classes
[~,Kp]=size(Ip) ; %
% Compute confusion matrix
C=zeros(M,M) ;    % Confusion matrix
for i=1:Kp
  m=Ix(Ip(1,i)) ; % true class in X
  j=Ip(2,i) ;     % classified class in Xp for same sample in X
  C(m,j)=C(m,j)+1 ;
end
Prei=zeros(1,M) ; % precision and recall for M classes
Reci=zeros(1,M) ;
for m=1:M
  TP=C(m,m) ;     % true positive
  FP=sum(C(:,m))-C(m,m) ;  % false positive
  FN=sum(C(m,:))-C(m,m) ;  % false negative
  Prei(m)=TP/(TP+FP) ;
  Reci(m)=TP/(TP+FN) ;
end
Pre=sum(Prei)/M ;
Rec=sum(Reci)/M ;
F1=2*Pre*Rec/(Pre+Rec) ;
% set NaN to 0
if isnan(Pre), Pre=0 ; end
if isnan(Rec), Rec=0 ; end
if isnan(F1), F1=0 ; end

end
