%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thumbnail: Fig. 4 without text and scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Intialization 
rng('default') ; 
key=1 ;               % key for random number generators
rng(key) ; 
%fprintf(1,'Emitter distance: %d (nm) \n',eD) ; 
%% Optical system 
na=1.40 ; 
lambda=723 ;            % Alexa700 wavelength in nm
alpha=2*pi*na/lambda ;  % =1.2167e-2
% 2D Gaussian PSF; sigma is estimated from Airy PSF
sigma=1.3238/alpha ;    % sigma=108.81; 2*sigma=217.61 (nm) 
FWHM=2*sqrt(2*log(2))*sigma ; % FWHM=256.22 (nm)
%% Frame 
% Region of view: [0,Lx]x[0,Ly]
Lx=2^12 ;
Ly=Lx ;               % frame size in nm
Dx=2^7 ; Dy=2^7 ;     % pixel size of cammera
Kx=Lx/Dx ; Ky=Ly/Dy ; % frame size in pixels
%% Emitter intensity and signal to noise ratio
Dt=0.01 ;             % s, time per frame (1/Dt is frame rate) 
Ih=300000 ;           % photons/s, average number of detected photons per second per emitter 
DtIh=Dt*Ih ;          % photon count per frame per emitter 
% 'mediumSNR'         % 
b=5 ;                 % photons/s/nm^2, Poisson noise
G=3 ;                 % photons/s/nm^2, variance of Gaussian noise
mu=0 ;                % photons/s/nm^2, mean of Gaussian noise 
rp=Ih/b ;             % 60000
rg=Ih/G ;             % 100000
r=rp*rg/(rp+rg) ;     % 37500
SNRdB=10*log10(r)-20*log10(sigma)-11.02 ; % -6.0127 (dB)
%% Emitter activations
M=150 ;               % # of emitters
N=200 ;               % Temporal resolution (TR): N*Dt=5 sec
K_ai=12 ;             % Average # of activations per emitter in data movie
                      % =(1-p0)*N ; ensure each emitter is activated at least once 
J=4 ;                 % Maximum state
r01=0.5 ;   r02=0.7 ;   r03=0.8 ;   r04=1.0 ; 
r21=1-r01 ; r32=1-r02 ; r43=1-r03 ;  
r00=1-K_ai/((N-K_ai)*(1+r21+r21*r32+r21*r32*r43)) ; % =0.9620
r10=1-r00 ; 
R=[r00 r01 r02 r03 r04    % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;                % =0.9400, probability of deactivation, i.e. state 0
p1=r10/den ;              % =0.0357, probability of state 1
p2=r10*r21/den ;          % =0.0179, probability of state 2
p3=r10*r21*r32/den ;      % =0.0054, probability of state 3
p4=r10*r21*r32*r43/den ;  % =0.0011, probability of state 4
pa=1-p0 ;                 % =0.0600, probability of activation 
Kaae=(1-p0)*M ;           % =K_a*M/K=9, average # of activated emitters/frame
pd=1-(1-p0^N)^M ;         % =6.3318e-04, probability that at least one emitter 
                          % is not activated in data movie
c0=zeros(M,N+1) ;         % states of Markov chains in data movie
for n=2:N+1
  for m=1:M
    c0(m,n)=(c0(m,n-1)+1)*(1-(rand<R(1,c0(m,n-1)+1))) ; % state transitions
  end
end
c=c0(:,2:N+1) ;           % remove initial all-0 state
if sum(c(:,1))==0         % avoid all zeros in 1st frame!
  c(1,1)=1 ; 
end
a=(c~=0) ;                % a(m,n)=1 if activated; a(m,n)=0 otherwise  
                          % sum(sum(a')==0): # of emitters never activated 
Kei=sum(a,1) ;            % number of activated emitters in nth frame
Kai=sum(a,2) ;            % number of activations for an emitter in movie
Ka=sum(Kai) ;             % total number of activations for all emitters in movie

%% Emitter locations - ground truth
D=[55 40 25] ;            % nm, distance between adjacent emitters 
h=zeros(size(D)) ;        % RMSE in theory
RMSE_=zeros(size(D)) ;    % RMSE by sample
RMSMD_=zeros(size(D)) ;   % RMSMD
RMSEP_=zeros(size(D)) ;   % RMSE_P
RMSMDP_=zeros(size(D)) ;  % RMSMD_P
% metrics after average of estimates for each emitter
h_a=zeros(size(D)) ;      % RMSE in theory
RMSE_a=zeros(size(D)) ;   % RMSE by sample
RMSMD_a=zeros(size(D)) ;  % RMSMD
RMSEP_a=zeros(size(D)) ;  % RMSE_P
RMSMDP_a=zeros(size(D)) ; % RMSMD_P
figure('Units','inches','Position',[2 2 3*2.3+0.04 3*2.3+0.04],'Color',[1 1 1]) ;
%% set image size and position
wx=2.48 ; wy=2.48 ; wxm=-0.15 ; wym=-0.1 ; 
x01=-0.14 ; x02=x01+wx-0.06 ; x03=x01+2*(wx-0.05)-0.18 ; 
y03=-0.125 ; y02=y03+wy-0.095 ; y01=y03+2*(wy-0.05)-0.20 ; 
for i=1:length(D)
  t1=(0:1/M:1-1/M) ;
  B=1500*t1+250 ;         % nm, distance from emitters to origin
  xy1=zeros(2,M) ;        % 2D emitter location
  xy1(:,1)=[0 ; B(1)] ;   % 1st emitter location
  for j=2:M
    theta0=atan2(xy1(2,j-1),xy1(1,j-1)) ; % angle of previous location
    d0=sqrt(xy1(:,j-1)'*xy1(:,j-1)) ;     % length of previous location
    d1=B(j) ; % distance of current location
    Dtheta=acos((d0^2+d1^2-D(i)^2)/(2*d0*d1)) ; % angle increment from previous to current location
    theta1=theta0-Dtheta ;  % angle of current location
    xy1(:,j)=[d1*cos(theta1) ; d1*sin(theta1)] ;
  end
  S=xy1+[0.5*Lx ; 0.5*Ly]*ones(1,M) ; % adjust to frame center
  Sn=zeros(2,max(Kei),N) ;  % activated emitter locations in nth frame
  pn=zeros(1,N) ;
  for n=1:N
    if Kei(n)>0
      for m=1:M
        if a(m,n)
          pn(n)=pn(n)+1 ;
          Sn(:,pn(n),n)=S(:,m) ;
        end
      end
    end
  end
  
  %% UGIA-F estimator 
  Xm=zeros(2,max(Kai),M) ;    % estimated locations for mth emitter
  CRLB=zeros(2,M) ;           % accumulative CRLBs for mth emitter
  pm=zeros(1,M) ;             % number of activations for mth emitter up to nth frame
  for n=1:N
    if Kei(n)>0
      [xyF,~,F_]=Gauss2D_UGIA_F(sigma,Kx,Ky,Dx,Dy,Dt,Ih,b,G,Sn(:,1:Kei(n),n)) ;
      h(i)=h(i)+sum(diag(F_)) ;
      p=0 ;
      for m=1:M
        if a(m,n)
          pm(m)=pm(m)+1 ; p=p+1 ;
          Xm(:,pm(m),m)=xyF(:,p) ;
          CRLB(1,m)=CRLB(1,m)+F_(2*(p-1)+1,2*(p-1)+1) ;
          CRLB(2,m)=CRLB(2,m)+F_(2*(p-1)+2,2*(p-1)+2) ;
        end
      end
    end
  end
  %% Compute metrics for UGIA-F image 
  X=zeros(2,Ka) ;             % all estimated locations from 1st to Mth emitters
  Ix=zeros(1,Ka) ;            % ith location in X is an estimate of Ix(i)th emitter
  p=0 ;
  for m=1:M
    X(:,p+1:p+Kai(m))=Xm(:,1:Kai(m),m) ;  % No action if Kai(m)=0
    Ix(p+1:p+Kai(m))=m ; 
    p=p+Kai(m) ;              % # of locations for emitters 1 up to m
  end
  % RMSE in theory
  h(i)=sqrt(h(i)/Ka) ; 
  % RMSE 
  RMSE_(i)=RMSE_P(S,X,Kai) ; 
  % RMSMD
  [RMSMD_(i),~]=RMSMD(S,X) ; 
  % partition
  [Xp,Kpi,Ip]=partitionX(S,X,Kai) ; 
  % F1 score
  [F1,Pre,Rec]=MCmetrics(Ip,Ix) ;
  % RMSE_P
  RMSEP_(i)=RMSE_P(S,Xp,Kpi) ; 
  % RMSMD_P
  RMSMDP_(i)=RMSMD_P(S,Xp,Kpi) ; 
  
  %% Compute metrics: after average over Xi 
  Ixa=zeros(1,M) ;            % ith location in Xa is an estimate of Ixa(i)th emitter
  for m=1:M
    h_a(i)=h_a(i)+sum(CRLB(:,m))/Kai(m)^2 ; % Note dividing by Kai(m) again
  end
  h_a(i)=sqrt(h_a(i)/M) ; 
  X_a=zeros(2,M) ;            % all estimated locations for M emitters
  for m=1:M
    for j=1:Kai(m)
      X_a(:,m)=X_a(:,m)+Xm(:,j,m) ;
    end
    X_a(:,m)=X_a(:,m)/Kai(m) ;  % averaging all estimates for an emitter 
    Ixa(m)=m ; 
  end
  Kaia=ones(1,M) ;
  % RMSE 
  RMSE_a(i)=RMSE_P(S,X_a,Kaia) ; 
  % RMSMD
  [RMSMD_a(i),~]=RMSMD(S,X_a) ; 
  % partition
  [Xp,Kpi,Ipa]=partitionX(S,X_a,Kaia) ; 
  % F1 score
  [F1a,Prea,Reca]=MCmetrics(Ipa,Ixa) ;
  % RMSE_P 
  RMSEP_a(i)=RMSE_P(S,Xp,Kpi) ; 
  % RMSMD_P
  RMSMDP_a(i)=RMSMD_P(S,Xp,Kpi) ; 
  fprintf(1,'D=%2d  F1=%6.3f h =%6.2f RMSE =%6.2f RMSEP =%6.2f RMSMDP =%6.2f RMSMD =%6.2f  \n', ...
    D(i),F1,h(i),RMSE_(i),RMSEP_(i),RMSMDP_(i),RMSMD_(i)) ;
  fprintf(1,'D=%2d  F1a=%5.3f ha=%6.2f RMSEa=%6.2f RMSEPa=%6.2f RMSMDPa=%6.2f RMSMDa=%6.2f  \n', ...
    D(i),F1a,h_a(i),RMSE_a(i),RMSEP_a(i),RMSMDP_a(i),RMSMD_a(i)) ;

  %% show images 
  switch i
    case 1
      % show 10th frame
      n=10 ;
      Fig1=subplot(3,3,1) ;
      U=Gauss2D_Frame(sigma,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,Sn(:,1:Kei(n),n)) ;
      show8bimage(U,'Yes','gray','No') ; hold on
      plot(S(1,:)/Dx+0.5,S(2,:)/Dy+0.5,'w.','MarkerSize',4) ;
      plot(Sn(1,1:Kei(n),n)/Dx+0.5,Sn(2,1:Kei(n),n)/Dy+0.5,'r.') ;
%     plot([200 700]/Dx+0.5,(Ly-[200 200])/Dy+0.5,'w-', ...  % scale bar = 500 nm
%        [200 200]/Dx+0.5,(Ly-[200-60 200+60])/Dy+0.5,'w-',[700 700]/Dx+0.5,(Ly-[200-60 200+60])/Dy+0.5,'w-') ; 
      hold off
%     text(200/Dx,(Ly-350)/Dy,'500 nm','Color','white','FontSize',8)
%     text(28.5,3,'(a1)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig1,'Units','inches','OuterPosition',[x01,y01,wx,wy]);
      % show X
      Fig4=subplot(3,3,4) ;     % show all estimated locations
      U=zeros(size(U)) ;
      show8bimage(U,'No','gray','No') ; hold on
      plot(X(1,:)/Dx+0.5,X(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
%     text(28.5,3,'(a2)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig4,'Units','inches','OuterPosition',[x01,y02,wx,wy+wym]);
      % show X_a
      Fig7=subplot(3,3,7) ;     % show all estimated locations
      show8bimage(U,'No','gray','No') ; hold on
      plot(X_a(1,:)/Dx+0.5,X_a(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
%     text(28.5,3,'(a3)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig7,'Units','inches','OuterPosition',[x01,y03,wx,wy]);
    case 2
      % show 10th frame
      n=10 ;
      Fig2=subplot(3,3,2) ;
      U=Gauss2D_Frame(sigma,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,Sn(:,1:Kei(n),n)) ;
      show8bimage(U,'Yes','gray','No') ; hold on
      plot(S(1,:)/Dx+0.5,S(2,:)/Dy+0.5,'w.','MarkerSize',4) ;
      plot(Sn(1,1:Kei(n),n)/Dx+0.5,Sn(2,1:Kei(n),n)/Dy+0.5,'r.') ; hold off
%     text(28.5,3,'(b1)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig2,'Units','inches','OuterPosition',[x02,y01,wx+wxm,wy]);
      % show X
      Fig5=subplot(3,3,5) ;     % show all estimated locations
      U=zeros(size(U)) ;
      show8bimage(U,'No','gray','No') ; hold on
      plot(X(1,:)/Dx+0.5,X(2,:)/Dy+0.5,'w.','MarkerSize',4) ; hold off
%     text(28.5,3,'(b2)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig5,'Units','inches','OuterPosition',[x02,y02,wx+wxm,wy+wym]);
      % show X_a
      Fig8=subplot(3,3,8) ;     % show all estimated locations
      show8bimage(U,'No','gray','No') ; hold on
      plot(X_a(1,:)/Dx+0.5,X_a(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
%     text(28.5,3,'(b3)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig8,'Units','inches','OuterPosition',[x02,y03,wx+wxm,wy]);
    case 3
      % show 10th frame
      n=10 ;
      Fig3=subplot(3,3,3) ;
      U=Gauss2D_Frame(sigma,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,Sn(:,1:Kei(n),n)) ;
      show8bimage(U,'Yes','gray','No') ; hold on
      plot(S(1,:)/Dx+0.5,S(2,:)/Dy+0.5,'w.','MarkerSize',4) ;
      plot(Sn(1,1:Kei(n),n)/Dx+0.5,Sn(2,1:Kei(n),n)/Dy+0.5,'r.') ; hold off
%     text(28.5,3,'(c1)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig3,'Units','inches','OuterPosition',[x03,y01,wx,wy]);
      % show X
      Fig6=subplot(3,3,6) ;     % show all estimated locations
      U=zeros(size(U)) ;
      show8bimage(U,'No','gray','No') ; hold on
      plot(X(1,:)/Dx+0.5,X(2,:)/Dy+0.5,'w.','MarkerSize',4) ; hold off
%     text(28.5,3,'(c2)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig6,'Units','inches','OuterPosition',[x03,y02,wx,wy+wym]);
      % show X_a
      Fig9=subplot(3,3,9) ;     % show all estimated locations
      show8bimage(U,'No','gray','No') ; hold on
      plot(X_a(1,:)/Dx+0.5,X_a(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
%     text(28.5,3,'(c3)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig9,'Units','inches','OuterPosition',[x03,y03,wx,wy]);
    otherwise
  end
  getframe(gcf) ;
end
%print('Fig4.tif','-dtiffn')
