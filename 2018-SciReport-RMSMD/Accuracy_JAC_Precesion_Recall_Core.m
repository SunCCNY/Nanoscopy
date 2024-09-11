% [accuracy JAC precision recall]=Accuracy_JAC_Precesion_Recall_Core(S,X,FWHM)
%
% Calculate accuracy, JAC, precision, and recall 
%
% Reference 
% [1] Sage, D. and et al. Quantitative evaluation of software packages for single-molecule 
% localization microscopy. Nat. Methods 12, 717-724 (2015).
%
% Input:
%   S         - Ground-truth fluorophore locations
%               size(S)=[dim,NS], NS dim-dimensional vectors, each represents a point
%   X         - Estimated fluorophore locations
%               size(X)=[dim,NX], NX dim-dimensional vectors,
%   FWHM      - FWHM of PSF
%
% Output:
%   accuracy  - also called Root mean square error (RMSE)
%   JAC       - Jaccard index, JAC=TP/(TP+FP+FN)
%   precision - precision=TP/(TP+FP)
%   recall    - recall=TP/(TP+FN)
%
% 04/08/2017, 05/04/2018
% Yi Sun

function [accuracy JAC precision recall]=Accuracy_JAC_Precesion_Recall_Core(S,X,FWHM)

[dim NS]=size(S) ; % dim - dimension of S, NS - # of points in S
[din NX]=size(X) ; % din - dimension of X, NX - # of points in X
if dim~=din,
  fprintf(1,'S and X have different dimensions!\n')
  return ;
end
if NS<=0||NX<=0,
  fprintf(1,'S and X cannot be empty!\n') ;
  return ;
end
% For each x in X, search s in S with minimum distance
D_X2S=zeros(1,NX) ; % minimum square distance from each X(:,i) to S
Xp=zeros(1,NX) ; % Xp(i)=j if X(:,i) is closest to S(:,j) than other S(:,k)  
for i=1:NX,
  tmp=S-X(:,i)*ones(1,NS) ;
  d2=sum(tmp.^2) ;
  [D_X2S(i) Xp(i)]=min(d2) ;
end
TP=0 ;
FP=0 ;
accuracy=0 ;
FWHM2=FWHM^2 ;
Sp=zeros(1,NS) ; % Sp(i)=1 if S(:,i)'s TP region contains an X(:,j); =0 otherwise 
for i=1:NX,
  if D_X2S(i)<=FWHM2,
    TP=TP+1 ;
    accuracy=accuracy+D_X2S(i) ;
    Sp(Xp(i))=1 ; 
  else
    FP=FP+1 ;
  end
end
accuracy=sqrt(accuracy/TP) ;
FN=sum(Sp==0) ;
precision=TP/(TP+FP) ;
recall=TP/(TP+FN) ;
JAC=TP/(TP+FP+FN) ;

end