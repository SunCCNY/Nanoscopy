%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates the figures in paper "Partition of Estimated 
% Locations: Approach to Quality Metric for Stochastic Optical 
% Localization Nanoscopy" 
% 
% References
% [1] Y. Sun, "Root mean square minimum distance as a quality metric for
% stochastic optical localization nanoscopy images," Sci. Reports, vol. 8, 
% no. 1, pp. 17211, Nov. 2018.
% [2] Y. Sun, "Localization precision of stochastic optical localization 
% nanoscopy using single frames," J. Biomed. Optics, vol. 18, no. 11, pp. 
% 111418-14, Oct. 2013.
% 
% Yi Sun
% Electrical Engineering Department
% The City College of City University of New York
% E-mail: ysun@ccny.cuny.edu
% 09/05/2021, 02/17/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 2. 8 emitters on a circle: FFL algorithm  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
rng('default') ; 
key=1 ;             % key for random number generators
rng(key) ; 
% set parameters
M=8 ;               % number of emitters on a circle
u=160 ;            % emitter distance 
Lx=500 ; Ly=Lx ;    % FOV size 
% generate emitter locations 
S=zeros(2,M) ;      % set of emitter locations
[S(:,1:M-1),diam]=circle2Dchart(u,M-1) ; 
% M-1 emitter locations on circle of radius = diam/2 = 184.38 nm ; 
% Half distance of S(:,1:M-1) to S(:,M) = diam/4 = 92.19 nm 
% -> Voronoi boundary between them 
S(:,1:M-1)=S(:,1:M-1)+[Lx/2 ; Ly/2]*ones(1,M-1) ; 
S(:,M)=[Lx/2 ; Ly/2] ;  % Mth emitter is at circle center 
rad=diam/2 ;        % distance between Mth and other emitters
% generate activations by Markov chain 
K_ai=15 ;           % average number of activations per emitter in data movie
J=4 ;               % Maximum state
r00=0.99 ;          % (1-p0)*M=11.85
r01=0.3 ;   r02=0.5 ;   r03=0.7 ;   r04=1.0 ; 
r10=1-r00 ; r21=1-r01 ; r32=1-r02 ; r43=1-r03 ;  
R=[r00 r01 r02 r03 r04  % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;          % =0.9789, probability of de-activation
pa=1-p0 ;           % =0.0211, probability of activation for an emitter in a frame
N=fix(K_ai/pa) ;     % =711, total number of frames in data movie 
c=zeros(M,N+1) ;    % states of Markov chains in data movie
c(1,1)=1 ;          % set 1st emitter activated in 1st frame
for n=2:N+1
  for m=1:M
    c(m,n)=(c(m,n-1)+1)*(1-(rand<R(1,c(m,n-1)+1))) ; % state transitions
  end
end
a=(c~=0) ;          % a(m,n)=1 if activated; a(m,n)=0 de-activated  
%% Frame-by-frame (FFL) localization 
Kai=sum(a,2) ;      % =17, 7, 21, 8, 11, 16, 12, 15, number of activations for ith emitter
Ka=sum(Kai) ;       % =107, total number of all activations in data movie
betai=Kai/sum(Kai) ;% = 0.1589, 0.0654, 0.1963, 0.0748, 0.1028, 0.1495, 0.1121, 0.1402
pa_=Ka/M/N ;        % =0.0188, estimated pa
pm=0.1 ;            % probability of missed activations per emitter by algorithm
Ki=fix((1-pm)*Kai); % =15, 6, 18, 7, 9, 14, 10, 13, number of estimates for ith emitter 
%Ki=Kai ;           % This causes little difference in performance 
%% Joint movie localization (JML)
%Ki=1*ones(1,M) ; 
%Ki(1)=2 ; Ki(5)=2 ; 
%% generate emitter estimates 
sigma=25 ;          % SD of estimates for an emitter 
K=sum(Ki) ;         % =92 total number of estimates for all emitters 
Xt=zeros(2,max(Ki),M) ;  % set of estimated locations 
X=zeros(2,K) ;
Ix=zeros(1,K) ;     % ith location in X is an estimate of Ix(i)th emitter
Xp=zeros(2,K) ;     % partition of X
pb=0:0.01:1 ;       % percentage of bias relative to rad
B=pb*rad ;          % bias distance from an emitter to Mth emitter
Pc=zeros(size(B)) ; % probability crosssing boundary of Voronoi cells 
RMSMD_=zeros(size(B)) ; 
RMSE_=zeros(size(B)) ;
RMSEP_=zeros(size(B)) ;  % RMSE after partition
RMSMDP_=zeros(size(B)) ; % RMSMD_P after partition
F1=zeros(size(B)) ;   % F1 score 
Pre=zeros(size(B)) ;  % Precision
Rec=zeros(size(B)) ;  % Recall 
mk =['r+' ; 'g+' ; 'b+' ; 'y+' ; 'c+' ; 'm+' ; 'w+' ; 'k+']  ;
mkp=['r.' ; 'g.' ; 'b.' ; 'y.' ; 'c.' ; 'm.' ; 'w.' ; 'k.']  ;
mks=['ro' ; 'go' ; 'bo' ; 'yo' ; 'co' ; 'mo' ; 'wo' ; 'ko']  ;
figure('Units','inches','Position',[1.5 2 3.4 3.4],'Color',[1 1 1]) ;
wx=1.73 ; wy=1.74 ; wxm=0.06 ; wym=-0.04 ; 
x01=-0.02 ; x02=x01+wx-0.06 ; 
y02=0.02 ; y01=y02+wy-0.06 ; 
[V1,V2]=voronoi(S(1,:),S(2,:)) ;  % Voronoi cells
for i=1:length(B)
  p=0 ; 
  for m=1:M
    Sbias=S(:,m)+(B(i)/rad)*(S(:,M)-S(:,m)) ;
    Xt(:,1:Ki(m),m)=sigma*randn(2,Ki(m))+Sbias ;
    X(:,p+1:p+Ki(m))=Xt(:,1:Ki(m),m) ;  % these are locations for mth emitter
    Ix(p+1:p+Ki(m))=m ; 
    p=p+Ki(m) ;
  end
  % probability crossing boundary 
  Pc(i)=1-Qfunc((B(i)-0.5*rad)/sigma) ; 
  % RMSE 
  RMSE_(i)=RMSE_P(S,X,Ki) ;
  % RMSMD
  [RMSMD_(i),~]=RMSMD(S,X) ;
  % partition
  [Xp,Kpi,Ip]=partitionX(S,X,Kai) ;
  % F1 score
  [F1(i),Pre(i),Rec(i)]=MCmetrics(Ip,Ix) ;
  % RMSE_P
  RMSEP_(i)=RMSE_P(S,Xp,Kpi) ;
  % RMSMD_P
  RMSMDP_(i)=RMSMD_P(S,Xp,Kpi) ;
  % (a) B(1)=0 -> B=0
  % Pc=Phi((B-0.5*rad)/sigma)=Phi(-0.5*rad/sigma)=Phi(-3.6876)=1-Qfunc(-3.6876)=1.1319e-4
  % F1(1)=100%
  if i==1    
    Fig2a=subplot(2,2,1);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(a)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig2a,'Units','inches','OuterPosition',[x01,y01,wx,wy]);
    set(gcf,'InvertHardcopy','off')
  end
  % (b) B(24)=42.4077 -> B=rad/2-2*sigma=42.1906
  % Pc=Phi((B-0.5*rad)/sigma)=Phi(-2)=1-Qfunc(-2)=0.0228
  % F1(24)=100%
  if i==24  
    Fig2b=subplot(2,2,2);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(b)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig2b,'Units','inches','OuterPosition',[x02,y01,wx+wxm,wy]);
    set(gcf,'InvertHardcopy','off')
  end
  % (c) B(51)=92.1906 -> B=rad/2=92.1906
  % Pc=Phi((B-0.5*rad)/sigma)=Phi(0)=1-Qfunc(0)=0.5
  % F1(51)=86.6%
  if i==51 
    Fig2c=subplot(2,2,3);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(c)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig2c,'Units','inches','OuterPosition',[x01,y02,wx,wy+wym]);
    set(gcf,'InvertHardcopy','off')
  end
  % (d) B(88)=160.4116 -> B=180 
  % Pc=Phi((180-0.5*rad)/sigma)=Phi(2.7124)=1-Qfunc(2.7124)=0.9967
  % F1(88)=31.3%
  if i==88 
    Fig2d=subplot(2,2,4);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(d)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig2d,'Units','inches','OuterPosition',[x02,y02,wx+wxm,wy+wym]);
    set(gcf,'InvertHardcopy','off')
  end
end
%print('Fig2a-d.tif','-dtiffn')

figure('Units','inches','Position',[3 4 3.4 2.5],'Color',[1 1 1]) ;
Fig2e=subplot(1,1,1);
h=(2*sigma^2+7*B.^2/8).^0.5 ;     % RMSE in theory
lg=plot(B,h,'m- .',B,RMSE_,'k- .',B,RMSEP_,'g- .',B,RMSMDP_,'r- .',B,RMSMD_,'b- .','MarkerSize',3.5) ; hold on
plot(B,100*Pc,'k--','MarkerSize',3.5) ; 
plot(B,100*F1,'k-','MarkerSize',3.5) ; 
% Note: Due to property of partition algorihtm, Prei=Reci in MCmetricsssss
% plot(B,100*Pre,'g-','MarkerSize',3.5) ;   
% plot(B,100*Rec,'b-','MarkerSize',3.5) ; 
plot(0*[1 1],[0 180],'k-.',(rad/2-2*sigma)*[1 1],[0 180],'k-.','MarkerSize',3.5) ; 
plot(rad/2*[1 1],[0 180],'k-.',160*[1 1],[0 180],'k-.','MarkerSize',3.5) ; hold off
axis([-1.5 161.5 0 160]) ;
grid on
legend(lg,'RMSE in theory, Eq. (17)','RMSE known {\itX_i}, Eq. (3)', ...
  'RMSE-P, Eq. (3)','RMSMD-P, Eq. (11)','RMSMD, Eq. (4)','Position',[0.27 0.722 0.1 0.2],'FontSize',5.5) ; 
text(142,105,'{\itP_c} (%)','FontSize',5.5) ; 
text(4,105,'F1 (%)','FontSize',5.5) ; 
text(2-0.5,10,'(a)','FontSize',5.5) ; 
text(45-0.5,10,'(b)','FontSize',5.5) ; 
text(95-0.5,10,'(c)','FontSize',5.5) ; 
text(153,10,'(d)','FontSize',5.5) ; 
set(gca,'XTick',[0 20 40 60 80 100 120 140 160 180]) ;
set(gca,'XTickLabel',[0 20 40 60 80 100 120 140 160 180]) ;
set(gca,'YTick',[0 20 40 60 80 100 120 140 160 180]) ;
set(gca,'YTickLabel',[0 20 40 60 80 100 120 140 160 180]) ;
text(142,155,'(e)','Color','black','FontSize',5.5)
xlabel('{\itB} (nm)','FontSize',5.5)
ylabel('Metrics (nm)','FontSize',5.5) 
set(gca,'FontSize',5.5) ;
set(Fig2e,'Units','inches','OuterPosition',[-0.1,0,3.7,2.6]);
set(gcf,'InvertHardcopy','off')
%print('Fig2e.tif','-dtiffn')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 3. 8 emitters on a circle: K~=M, JML algorithm  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
rng('default') ; 
key=1 ;             % key for random number generators
rng(key) ; 
% set parameters
M=8 ;               % number of emitters on a circle
u=160 ;             % emitter distance 
Lx=500 ; Ly=Lx ;    % FOV size 
% generate emitter locations 
S=zeros(2,M) ;      % set of emitter locations
[S(:,1:M-1),diam]=circle2Dchart(u,M-1) ; 
% M-1 emitter locations on circle of radius = diam/2 = 184.38 nm ; 
% Half distance of S(:,1:M-1) to S(:,M) = diam/4 = 92.19 nm 
% -> Voronoi boundary between them 
S(:,1:M-1)=S(:,1:M-1)+[Lx/2 ; Ly/2]*ones(1,M-1) ; 
S(:,M)=[Lx/2 ; Ly/2] ;  % Mth emitter is at circle center 
rad=diam/2 ;        % distance between Mth and other emitters
% generate activations by Markov chain 
K_ai=15 ;           % average number of activations per emitter in data movie
J=4 ;               % Maximum state
r00=0.99 ;          % (1-p0)*M=11.85
r01=0.3 ;   r02=0.5 ;   r03=0.7 ;   r04=1.0 ; 
r10=1-r00 ; r21=1-r01 ; r32=1-r02 ; r43=1-r03 ;  
R=[r00 r01 r02 r03 r04  % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;          % =0.9789, probability of de-activation
pa=1-p0 ;           % =0.0211, probability of activation for an emitter in a frame
N=fix(K_ai/pa) ;    % =711, total number of frames in data movie 
c=zeros(M,N+1) ;    % states of Markov chains in data movie
c(1,1)=1 ;          % set 1st emitter activated in 1st frame
for n=2:N+1
  for m=1:M
    c(m,n)=(c(m,n-1)+1)*(1-(rand<R(1,c(m,n-1)+1))) ; % state transitions
  end
end
a=(c~=0) ;          % a(m,n)=1 if activated; a(m,n)=0 de-activated  
%% Frame-by-frame (FFL) localization 
Kai=sum(a,2) ;      % =17, 7, 21, 8, 11, 16, 12, 15, number of activations for ith emitter
Ka=sum(Kai) ;       % =107, total number of all activations in data movie
%betai=Kai/sum(Kai) ;% = 0.1589, 0.0654, 0.1963, 0.0748, 0.1028, 0.1495, 0.1121, 0.1402
%pa_=Ka/M/N ;       % =0.0188, estimated pa
%pm=0.1 ;           % probability of missed activations per emitter by algorithm
%Ki=fix((1-pm)*Kai);% =15, 6, 18, 7, 9, 14, 10, 13, number of estimates for ith emitter 
%Ki=Kai ;           % This causes little difference in performance 
%% Joint movie localization (JML)
Ki=1*ones(1,M) ; 
Ki(1)=2 ; Ki(5)=2 ; 
%% generate emitter estimates 
sigma=25 ;          % SD of estimates for an emitter 
K=sum(Ki) ;         % =92 total number of estimates for all emitters 
Xt=zeros(2,max(Ki),M) ;  % set of estimated locations 
X=zeros(2,K) ;
Ix=zeros(1,K) ;     % ith location in X is an estimate of Ix(i)th emitter
Xp=zeros(2,K) ;     % partition of X
pb=0:0.01:1 ;       % percentage of bias relative to rad
B=pb*rad ;          % bias distance from an emitter to Mth emitter
Pc=zeros(size(B)) ; % probability crosssing boundary of Voronoi cells 
RMSMD_=zeros(size(B)) ;
RMSE_=zeros(size(B)) ;
RMSEP_=zeros(size(B)) ;  % RMSE after partition
RMSMDP_=zeros(size(B)) ; % RMSMD_P after partition
F1=zeros(size(B)) ;   % F1 score 
Pre=zeros(size(B)) ;  % Precision
Rec=zeros(size(B)) ;  % Recall 
mk =['r+' ; 'g+' ; 'b+' ; 'y+' ; 'c+' ; 'm+' ; 'w+' ; 'k+']  ;
mkp=['r.' ; 'g.' ; 'b.' ; 'y.' ; 'c.' ; 'm.' ; 'w.' ; 'k.']  ;
mks=['ro' ; 'go' ; 'bo' ; 'yo' ; 'co' ; 'mo' ; 'wo' ; 'ko']  ;
figure('Units','inches','Position',[1.5 2 3.4 3.4],'Color',[1 1 1]) ;
wx=1.73 ; wy=1.74 ; wxm=0.06 ; wym=-0.04 ; 
x01=-0.02 ; x02=x01+wx-0.06 ; 
y02=0.02 ; y01=y02+wy-0.06 ; 
[V1,V2]=voronoi(S(1,:),S(2,:)) ;  % Voronoi cells
for i=1:length(B)
  p=0 ;
  for m=1:M
    Sbias=S(:,m)+(B(i)/rad)*(S(:,M)-S(:,m)) ;
    Xt(:,1:Ki(m),m)=sigma*randn(2,Ki(m))+Sbias ;
    X(:,p+1:p+Ki(m))=Xt(:,1:Ki(m),m) ;
    Ix(p+1:p+Ki(m))=m ; 
    p=p+Ki(m) ;
  end
  % probability crossing boundary 
  Pc(i)=1-Qfunc((B(i)-0.5*rad)/sigma) ; 
  % RMSE
  RMSE_(i)=RMSE_P(S,X,Ki) ;
  % RMSMD
  [RMSMD_(i),~]=RMSMD(S,X) ;
  % partition
  [Xp,Kpi,Ip]=partitionX(S,X,Kai) ;
  % F1 score
  [F1(i),Pre(i),Rec(i)]=MCmetrics(Ip,Ix) ;
  % RMSE_P
  RMSEP_(i)=RMSE_P(S,Xp,Kpi) ;
  % RMSMD_P
  RMSMDP_(i)=RMSMD_P(S,Xp,Kpi) ;
  % (a) D(1)=0 -> B=0
  % Pc=Phi((B-0.5*rad)/sigma)=Phi(-0.5*rad/sigma)=Phi(-3.6876)=1-Qfunc(-3.6876)=1.1319e-4
  % F1(1)=100%
  if i==1    
    Fig3a=subplot(2,2,1);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(a)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig3a,'Units','inches','OuterPosition',[x01,y01,wx,wy]);
    set(gcf,'InvertHardcopy','off')
  end
  % (b) D(24)=42.4077 -> B=rad/2-2*sigma=42.1906
  % Pc=Phi((B-0.5*rad)/sigma)=Phi(-2)=1-Qfunc(-2)=0.0228
  % F1(24)=100%
  if i==24  
    Fig3b=subplot(2,2,2);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(b)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig3b,'Units','inches','OuterPosition',[x02,y01,wx+wxm,wy]);
    set(gcf,'InvertHardcopy','off')
  end
  % (c) D(51)=92.1906 -> B=rad/2=92.1906
  % Pc=Phi((B-0.5*rad)/sigma)=Phi(0)=1-Qfunc(0)=0.5
  % F1(51)=93.8%
  if i==51 
    Fig3c=subplot(2,2,3);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(c)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig3c,'Units','inches','OuterPosition',[x01,y02,wx,wy+wym]);
    set(gcf,'InvertHardcopy','off')
  end
  % (d) D(88)=160.4116 -> B=160 
  % Pc=Phi((160-0.5*rad)/sigma)=Phi(2.7124)=1-Qfunc(2.7124)=0.9967
  % F1(88)=62.5%
  if i==88 
    Fig3d=subplot(2,2,4);
    plot(1e-3*V1,1e-3*V2,'k:','MarkerSize',3.5)  ; hold on
    % show: estimates with bias
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*X(1,p+1:p+Ki(m)),1e-3*X(2,p+1:p+Ki(m)),mk(j,:),'MarkerSize',3.5) ;
      p=p+Ki(m) ;
    end
    p=0 ; 
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*Xp(1,p+1:p+Kpi(m)),1e-3*Xp(2,p+1:p+Kpi(m)),mkp(j,:),'MarkerSize',3.5) ;
      p=p+Kpi(m) ;
    end
    % show: emitter locations
    for m=1:M
      j=mod(m-1,8)+1 ;
      plot(1e-3*S(1,m),1e-3*S(2,m),mks(j,:),'MarkerSize',3.5) ;
    end
    hold off
    axis(1e-3*[0 Lx 0 Ly]) ;
    set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'XTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]) ;
    set(gca,'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5]) ;
    text(0.42,0.475,'(d)','Color','black','FontSize',5.5)
    xlabel('{\mu}m','FontSize',5.5) ;
    ylabel('{\mu}m','FontSize',5.5) ;
    set(gca,'color',[0.75 0.75 0.75])
    set(gca,'FontSize',5.5) ; 
    set(Fig3d,'Units','inches','OuterPosition',[x02,y02,wx+wxm,wy+wym]);
    set(gcf,'InvertHardcopy','off')
  end
end
%print('Fig3a-d.tif','-dtiffn')

figure('Units','inches','Position',[3 4 3.4 2.5],'Color',[1 1 1]) ;
Fig3e=subplot(1,1,1);
h=(2*sigma^2+7*B.^2/8).^0.5 ;       % RMSE in theory
lg=plot(B,h,'m- .',B,RMSE_,'k- .',B,RMSEP_,'g- .',B,RMSMDP_,'r- .',B,RMSMD_,'b- .','MarkerSize',3.5) ; hold on
plot(B,100*Pc,'k--','MarkerSize',3.5) ; 
plot(B,100*F1,'k-','MarkerSize',3.5) ; 
% Note: Due to property of partition algorihtm, Reci=Nan in MCmetricsssss
% plot(B,100*Pre,'g-','MarkerSize',3.5) ;   
% plot(B,100*Rec,'b-','MarkerSize',3.5) ; 
plot(0*[1 1],[0 180],'k-.',(rad/2-2*sigma)*[1 1],[0 180],'k-.','MarkerSize',3.5) ; 
plot(rad/2*[1 1],[0 180],'k-.',160*[1 1],[0 180],'k-.','MarkerSize',3.5) ; hold off
axis([-1.5 161.5 0 160]) ;
grid on
legend(lg,'RMSE in theory, Eq. (17)','RMSE known {\itX_i}, Eq. (3)', ...
  'RMSE-P, Eq. (3)','RMSMD-P, Eq. (11)','RMSMD, Eq. (4)','Position',[0.27 0.722 0.1 0.2],'FontSize',5.5) ; 
text(142,105,'{\itP_c} (%)','FontSize',5.5) ; 
text(4,105,'F1 (%)','FontSize',5.5) ; 
text(2-0.5,10,'(a)','FontSize',5.5) ; 
text(45-0.5,10,'(b)','FontSize',5.5) ; 
text(95-0.5,10,'(c)','FontSize',5.5) ; 
text(153,10,'(d)','FontSize',5.5) ; 
set(gca,'XTick',[0 20 40 60 80 100 120 140 160 180]) ;
set(gca,'XTickLabel',[0 20 40 60 80 100 120 140 160 180]) ;
set(gca,'YTick',[0 20 40 60 80 100 120 140 160 180]) ;
set(gca,'YTickLabel',[0 20 40 60 80 100 120 140 160 180]) ;
text(142,155,'(e)','Color','black','FontSize',5.5)
xlabel('{\itB} (nm)','FontSize',5.5)
ylabel('Metrics (nm)','FontSize',5.5) 
set(gca,'FontSize',5.5) ;
set(Fig3e,'Units','inches','OuterPosition',[-0.1,0,3.7,2.6]);
set(gcf,'InvertHardcopy','off')
%print('Fig3e.tif','-dtiffn')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 4 UGIA-F estimator  
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
      plot([200 700]/Dx+0.5,(Ly-[200 200])/Dy+0.5,'w-', ...  % scale bar = 500 nm
        [200 200]/Dx+0.5,(Ly-[200-60 200+60])/Dy+0.5,'w-',[700 700]/Dx+0.5,(Ly-[200-60 200+60])/Dy+0.5,'w-') ; 
      hold off
      text(200/Dx,(Ly-350)/Dy,'500 nm','Color','white','FontSize',8)
      text(28.5,3,'(a1)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig1,'Units','inches','OuterPosition',[x01,y01,wx,wy]);
      % show X
      Fig4=subplot(3,3,4) ;     % show all estimated locations
      U=zeros(size(U)) ;
      show8bimage(U,'No','gray','No') ; hold on
      plot(X(1,:)/Dx+0.5,X(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
      text(28.5,3,'(a2)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig4,'Units','inches','OuterPosition',[x01,y02,wx,wy+wym]);
      % show X_a
      Fig7=subplot(3,3,7) ;     % show all estimated locations
      show8bimage(U,'No','gray','No') ; hold on
      plot(X_a(1,:)/Dx+0.5,X_a(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
      text(28.5,3,'(a3)','Color','white','fontsize',8) ;
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
      text(28.5,3,'(b1)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig2,'Units','inches','OuterPosition',[x02,y01,wx+wxm,wy]);
      % show X
      Fig5=subplot(3,3,5) ;     % show all estimated locations
      U=zeros(size(U)) ;
      show8bimage(U,'No','gray','No') ; hold on
      plot(X(1,:)/Dx+0.5,X(2,:)/Dy+0.5,'w.','MarkerSize',4) ; hold off
      text(28.5,3,'(b2)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig5,'Units','inches','OuterPosition',[x02,y02,wx+wxm,wy+wym]);
      % show X_a
      Fig8=subplot(3,3,8) ;     % show all estimated locations
      show8bimage(U,'No','gray','No') ; hold on
      plot(X_a(1,:)/Dx+0.5,X_a(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
      text(28.5,3,'(b3)','Color','white','fontsize',8) ;
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
      text(28.5,3,'(c1)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig3,'Units','inches','OuterPosition',[x03,y01,wx,wy]);
      % show X
      Fig6=subplot(3,3,6) ;     % show all estimated locations
      U=zeros(size(U)) ;
      show8bimage(U,'No','gray','No') ; hold on
      plot(X(1,:)/Dx+0.5,X(2,:)/Dy+0.5,'w.','MarkerSize',4) ; hold off
      text(28.5,3,'(c2)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig6,'Units','inches','OuterPosition',[x03,y02,wx,wy+wym]);
      % show X_a
      Fig9=subplot(3,3,9) ;     % show all estimated locations
      show8bimage(U,'No','gray','No') ; hold on
      plot(X_a(1,:)/Dx+0.5,X_a(2,:)/Dy+0.5,'w.','MarkerSize',4) ;  hold off
      text(28.5,3,'(c3)','Color','white','fontsize',8) ;
      axis off
      set(gca,'XTick',[]) ;     % Turn off X and Y ticks
      set(gca,'YTick',[]) ;
      set(Fig9,'Units','inches','OuterPosition',[x03,y03,wx,wy]);
    otherwise
  end
  getframe(gcf) ;
end
%print('Fig4.tif','-dtiffn')

%% result: Na=1806 
% D=55  F1= 0.984 h = 16.54 RMSE = 16.17 RMSEP = 15.87 RMSMDP = 15.29 RMSMD = 14.74  
% D=55  F1a=1.000 ha=  5.40 RMSEa=  5.04 RMSEPa=  5.04 RMSMDPa=  5.04 RMSMDa=  5.04  
% D=40  F1= 0.914 h = 28.60 RMSE = 29.66 RMSEP = 29.39 RMSMDP = 28.26 RMSMD = 24.67  
% D=40  F1a=1.000 ha=  9.67 RMSEa= 10.20 RMSEPa= 10.20 RMSMDPa= 10.20 RMSMDa=  8.85  
% D=25  F1= 0.711 h =149.99 RMSE =207.95 RMSEP =202.13 RMSMDP =194.22 RMSMD =159.43  
% D=25  F1a=0.840 ha= 50.39 RMSEa= 72.27 RMSEPa= 70.39 RMSMDPa= 70.39 RMSMDa= 47.17 
