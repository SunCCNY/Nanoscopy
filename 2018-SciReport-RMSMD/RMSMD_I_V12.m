%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Root mean square mimimum distance (RMSMD) as a quality metric 
% for localization nanoscopy images
%
% Yi Sun
% Electrical Engineering Department
% The City College of City University of New York
% E-mail: ysun@ccny.cuny.edu
% 03/14/2017, 04/04/2017, 09/12/2017, 11/18/2017, 01/19/2018, 
% 05/08/2018, 06/03/2018, 08/21/2018, 09/11/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 1. Voronoi cells 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% S, Voronoi cells
dim=2 ; 
S=[200 750 1190 1200
   460 660  600  200] ; 
Fig1=figure('Position',[400 300 500-20 500*(1000/1400)*(102/108)],'Color',[1 1 1]) ;
plot(S(1,:),S(2,:),'ko','MarkerSize',4) ; hold on
axis([0 1400 0 1000]) ;
[V1 V2]=voronoi(S(1,:),S(2,:)) ;
plot(V1,V2,'k--')  ; 
% X, Voronoi cells
X=[300  260  550 820  950 1100
   400 760 800 560 800  100] ; 
plot(X(1,:),X(2,:),'k.') ; 
[V1 V2]=voronoi(X(1,:),X(2,:)) ;
plot(V1,V2,'k:')  ; 
% Denote S
text(120,460,'{\its}_1') ; 
text(690,600,'{\its}_2') ; 
text(1210,630,'{\its}_3') ; 
text(1240,240,'{\its}_4') ; 
% Denote X
text(280,340,'{\itx}_1') ; 
text(190,770,'{\itx}_2') ; 
text(530,860,'{\itx}_3') ; 
text(790,500,'{\itx}_4') ; 
text(910,860,'{\itx}_5') ; 
text(1050,50,'{\itx}_6') ; 
% Arrows
s=[220 ; 450] ; x=[285 ; 411] ;
plot([s(1) x(1)],[s(2) x(2)],'k-') ; 
draw2DArrow(s,x,30*ones(2,1),0.12*pi,'k-',0.5) ; 
draw2DArrow(x,s,30*ones(2,1),0.1*pi,'k-',0.5) ; 
s=[204 ; 490] ; x=[255 ; 740] ;
plot([s(1) x(1)],[s(2) x(2)],'k-') ; 
draw2DArrow(x,s,30*ones(2,1),0.1*pi,'k-',0.5) ; 
s=[728 ; 678] ; x=[570 ; 790] ;
plot([s(1) x(1)],[s(2) x(2)],'k-') ; 
draw2DArrow(x,s,30*ones(2,1),0.1*pi,'k-',0.5) ; 
s=[765 ; 637] ; x=[810 ; 575] ;
plot([s(1) x(1)],[s(2) x(2)],'k-') ; 
draw2DArrow(s,x,30*ones(2,1),0.1*pi,'k-',0.5) ; 
draw2DArrow(x,s,30*ones(2,1),0.1*pi,'k-',0.5) ; 
s=[1183 ; 183] ; x=[1115 ; 115] ;
plot([s(1) x(1)],[s(2) x(2)],'k-') ; 
draw2DArrow(s,x,30*ones(2,1),0.1*pi,'k-',0.5) ; 
draw2DArrow(x,s,30*ones(2,1),0.1*pi,'k-',0.5) ; 
s=[775 ; 678] ; x=[928 ; 785] ;
plot([s(1) x(1)],[s(2) x(2)],'k-') ; 
draw2DArrow(x,s,30*ones(2,1),0.1*pi,'k-',0.5) ; 
s=[1167 ; 622] ; x=[970 ; 790] ;
plot([s(1) x(1)],[s(2) x(2)],'k-') ; 
draw2DArrow(s,x,30*ones(2,1),0.1*pi,'k-',0.5) ; 
hold off
%grid on
xlabel('{\itx} (nm)') ; 
ylabel('{\ity} (nm)') ;
set(gca,'XTick',[0 200 400 600 800 1000 1200 1400]) ; 
set(gca,'XTickLabel',[0 200 400 600 800 1000 1200 1400]') ; 
set(gca,'YTick',[0 200 400 600 800 1000]) ; 
set(gca,'YTickLabel',[0 200 400 600 800 1000]') ; 
AxesHandle=findobj(Fig1,'Type','axes') ; 
set(AxesHandle,'Position',[0.115,0.11,0.85,0.85]);
%print -deps Fig1 % (to keep physical size, in Figure window, use File -> Save As)
 
[RMSMD0 tmp]=RMSMD(S,X) ; 
[RMSMD0 RMSMD0^2] % = 201.84 nm ; RMSMD^2 = 40740 nm^2
[RMSMDX_ tmp]=RMSMD(S(:,1:2),[X(:,2) X(:,3)]) ; 
[RMSMDX_ RMSMDX_^2] % = 276.76 nm ; RMSMDX_^2 = 76600 nm^2
[RMSMDX_X tmp]=RMSMD(S,[X(:,1) X(:,4:6)]) ; 
[RMSMDX_X RMSMDX_X^2] % = 178.26 nm ; RMSMDX_X^2 = 31775 nm^2
[RMSMDXkSk tmp]=RMSMD([S(:,1:2) S(:,4)],[X(:,1) X(:,4) X(:,6)]) ; 
[RMSMDXkSk RMSMDXkSk^2] % = 127.15 nm ; RMSMDXkSk^2 = 16166.67 nm^2

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 2. X0 is fixed and X1 changes on a trajectory
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
r=750 ; % radius of circle on which si's are located 
M=7 ;   % # of si
% S
theta=2*pi*(0:M-1)/M ;
S=r*[cos(theta) ; sin(theta)] ;
Fig2=figure('Position',[100 50 620 600],'Color',[1 1 1]) ;
% 2*r*sin(2*pi/(2*M)) = 650.82 nm -- distance between si's
% Voronoi cells of si's
Fig2a=subplot(2,2,1) ;
line=3*r*[cos(theta+2*pi/(2*M)) ; sin(theta+2*pi/(2*M))] ;
for i=1:M,
  plot([0 line(1,i)],[0,line(2,i)],'k--','LineWidth',0.125) ; hold on
end
% X consists of M fixed X0, and M X1 moving along their trajectories 
X0=1.35*S ; 
I=(0:200)/200 ; 
ii=20 ; 
alpha=0.8*2*pi*I ;
X10=1300*[I.*cos(alpha) ; I.*sin(alpha)] ; % trajectory for xi, theta=0
X10=X10+[2500*I ; zeros(size(I))] ; 
X1=zeros(2*M,length(I)) ; % X1 is changing
for i=1:M,
  T=[cos(theta(i)) -sin(theta(i))
     sin(theta(i))  cos(theta(i))] ;
  X1(2*(i-1)+1:2*(i-1)+2,:)=T*X10 ; 
  plot(X1(2*(i-1)+1,:),X1(2*(i-1)+2,:),'k-','LineWidth',0.125) ; 
  plot(X1(2*(i-1)+1,ii),X1(2*(i-1)+2,ii),'k.','MarkerSize',5) ; 
  plot(X0(1,i),X0(2,i),'k.','MarkerSize',5) ; 
end
% Draw arrows on X1
jj=155 ; 
a0=zeros(2,M) ; a1=zeros(2,M) ; 
for i=1:M, 
  a0(:,i)=X1(2*(i-1)+1:2*(i-1)+2,jj-2) ; 
  a1(:,i)=X1(2*(i-1)+1:2*(i-1)+2,jj) ; 
end
draw2DArrow(a0,a1,80*ones(2,1),0.15*pi,'k-',0.125) ; 
hold on
% Circle with radius of PSF FWHM
FWHM=184 ; % FWHM of PSF = 2*sqrt(2*log(2))*78.26 = 184.288
t=(0:100)/100 ; 
cc0=FWHM*[cos(2*pi*t) ; sin(2*pi*t)] ; 
for i=1:M,
  cc=S(:,i)*ones(1,length(t))+cc0 ; 
  plot(cc(1,:),cc(2,:),'k:','LineWidth',0.125) ; 
end
plot(S(1,:),S(2,:),'ko','MarkerSize',4) ; 
hold off
axis([-1500 1500 -1500 1500]) ;
text(-1400,1350,'(a)') ; 
xlabel('{\itx} (nm)','FontSize',8) ; 
ylabel('{\ity} (nm)','FontSize',8) ; 
set(gcf,'Color',[1 1 1])
% Denote S
text(600,-10,'{\its}_1','FontSize',8) ; 
text(390,490,'{\its}_2','FontSize',8) ; 
text(-185,635,'{\its}_3','FontSize',8) ; 
text(-630,315,'{\its}_4','FontSize',8) ; 
text(-640,-275,'{\its}_5','FontSize',8) ; 
text(-180,-630,'{\its}_6','FontSize',8) ; 
text(332,-510,'{\its}_7','FontSize',8) ; 
% Denote X0
text(1060,0,'{\itx}_1','FontSize',8) ; 
text(660,840,'{\itx}_2','FontSize',8) ; 
text(-290,1100,'{\itx}_3','FontSize',8) ; 
text(-1080,470,'{\itx}_4','FontSize',8) ; 
text(-1040,-470,'{\itx}_5','FontSize',8) ; 
text(-280,-1090,'{\itx}_6','FontSize',8) ; 
text(670,-855,'{\itx}_7','FontSize',8) ; 
% Denote X1
text(340,-35,'{\itx}_8','FontSize',8) ; 
text(220,300,'{\itx}_9','FontSize',8) ; 
text(-110,380,'{\itx}_{10}','FontSize',8) ; 
text(-370,195,'{\itx}_{11}','FontSize',8) ; 
text(-500,-120,'{\itx}_{12}','FontSize',8) ; 
text(-210,-355,'{\itx}_{13}','FontSize',8) ; 
text(175,-310,'{\itx}_{14}','FontSize',8) ; 
set(gca,'FontSize',8) ; 
set(Fig2a,'OuterPosition',[0,0.5,0.5,0.5]); 

%% calculate metrics
RMSMD0=zeros(size(I)) ; 
accuracy=zeros(size(I)) ;
JAC=zeros(size(I)) ; 
precision=zeros(size(I)) ; 
recall=zeros(size(I)) ; 
d=zeros(size(I)) ; 
XX=zeros(size(S)) ; 
for i=1:length(I),
  d(i)=sqrt(X1(1:2,i)'*X1(1:2,i)) ; 
  for j=1:M,
    XX(:,j)=X1(2*(j-1)+1:2*(j-1)+2,i) ; 
  end
  [RMSMD0(i) tmp]=RMSMD(S,[X0 XX]) ; 
  [accuracy(i) JAC(i) precision(i) recall(i)]=Accuracy_JAC_Precesion_Recall_Core(S,[X0 XX],FWHM) ;
end

Fig2b=subplot(2,2,2) ;
  plot(d,RMSMD0,'k-','LineWidth',0.125) ; hold on
  plot(d(ii),RMSMD0(ii),'k.','MarkerSize',5) ; 
  axis([0 1500 0 600]) ;
  text(60,570,'(b)') ; 
  xlabel('{\itd} (nm)','FontSize',8) ; 
  ylabel('RMSMD (nm)','FontSize',8) ; 
set(gca,'XTick',[0 250 500 750 1000 1250 1500]) ; 
set(gca,'XTickLabel',[0 250 500 750 1000 1250 1500]) ; 
set(gca,'YTick',[0 100 200 300 400 500 600]) ; 
set(gca,'YTickLabel',[0 100 200 300 400 500 600]) ; 
set(gca,'FontSize',8) ; 
set(Fig2b,'OuterPosition',[0.51,0.5,0.515,0.5]); 

Fig2c=subplot(2,2,3) ;
  plot(d,accuracy,'k-','LineWidth',0.125) ; hold on
  plot(d(ii),accuracy(ii),'k.','MarkerSize',5) ; hold off
  axis([0 1500 0 600]) ;
  text(60,570,'(c)') ; 
  xlabel('{\itd} (nm)','FontSize',8) ; 
  ylabel('Accuracy(nm)','FontSize',8) ; 
set(gca,'XTick',[0 250 500 750 1000 1250 1500]) ; 
set(gca,'XTickLabel',[0 250 500 750 1000 1250 1500]) ; 
set(gca,'YTick',[0 100 200 300 400 500 600]) ; 
set(gca,'YTickLabel',[0 100 200 300 400 500 600]) ; 
set(gca,'FontSize',8) ; 
set(Fig2c,'OuterPosition',[0.016,0,0.486,0.486]); 

Fig2d=subplot(2,2,4) ;
mk=[5 33 63 128 151] ; 
  plot(d,100*precision,'k-','LineWidth',0.125) ; hold on
	plot(d(mk),100*precision(mk),'ko','MarkerSize',4) ;   
  plot(d(ii),100*precision(ii),'k.','MarkerSize',5) ; 
lg(1)=plot([-20 -10],[-20 -10],'k-o','LineWidth',0.125,'MarkerSize',4) ;   
%
mk=[9 39 140 129 155] ; 
  plot(d,100*recall,'k-','LineWidth',0.125) ; hold on
	plot(d(mk),100*recall(mk),'ks','MarkerSize',4) ;   
  plot(d(ii),100*recall(ii),'k.','MarkerSize',5) ; 
lg(2)=plot([-20 -10],[-20 -10],'k-s','LineWidth',0.125,'MarkerSize',4) ;   
%
mk=[7 36 74 131 153] ; 
  plot(d,100*JAC,'k-','LineWidth',0.125) ; hold on
  plot(d(mk),100*JAC(mk),'kx','MarkerSize',6) ;   
  plot(d(ii),100*JAC(ii),'k.','MarkerSize',5) ; 
lg(3)=plot([-20 -10],[-20 -10],'k-x','LineWidth',0.125,'MarkerSize',6) ;   
  hold off
  axis([0 1500 -10 120]) ;
  legend(lg,'Precision','Recall','JAC','location','northeast') ;
  text(60,113,'(d)') ; 
  xlabel('{\itd} (nm)','FontSize',8) ; 
  ylabel('Precision, Recall, JAC (%)','FontSize',8) ; 
set(gca,'XTick',[0 250 500 750 1000 1250 1500]) ; 
set(gca,'XTickLabel',[0 250 500 750 1000 1250 1500]) ; 
set(gca,'YTick',[0 20 40 60 80 100 120]) ; 
set(gca,'YTickLabel',[0 20 40 60 80 100 120]) ; 
set(gca,'FontSize',8) ; 
set(Fig2d,'OuterPosition',[0.51,0,0.515,0.486]); 
%print -deps Fig2 % (to keep physical size, in Figure window, use File -> Save As)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 3-5. Helix
%  2-D emitters are located on a helix
clear
% optical system parameters
na=1.4 ; lambda=520 ; % nm
a=2*pi*na/lambda ;
% estimation of sigma of Airy PSF
rc=7.016/a ;
d=1 ; r=0.01:d:rc ; 
sigma=d*sqrt(8/pi)*sum(besselj(1,a*r).^2)  % = 78.26 nm
FWHM=2*sqrt(2*log(2))*sigma ; % = 184.28 nm
% emitter optics
Dt=0.01 ;     % second, time per frame (1/Dt is frame rate) 
Ih=300000 ;   % average number of detected photons per emitter per second
DtIh=Dt*Ih ;  % photon count per frame per emitter 
SPNR=0.2 ; SGNR=0.3 ; % um^2/emitter
SNR=SPNR*SGNR/(SPNR+SGNR) ;
% sample is located at [0,Lx]x[0,Ly]
Lx=2^11 ; Ly=2^11 ;   % frame size in nm
Dx=128 ; Dy=128 ;     % pixel size of cammera
Kx=Lx/Dx ; Ky=Ly/Dy ; % frame size in pixels
K=Kx*Ky ;							% total number of pixels per frame
% emitter locations
M=250 ; % # of emitters
Dd=25 ; % nm, distance between adjacent emitters 
t1=0:1/M:1-1/M ; 
D=(1000*t1).^0.97+40 ; % nm, distance from emitters to origin
xy=zeros(2,M) ; % 2D emitter locations
xy(:,1)=[0 ; D(1)] ; 
for i=2:M, 
  theta0=atan2(xy(2,i-1),xy(1,i-1)) ; % angle of previous location
  d0=sqrt(xy(:,i-1)'*xy(:,i-1)) ;     % length of previous location
  d1=D(i) ; % distance of current location
  Dtheta=acos((d0^2+d1^2-Dd^2)/(2*d0*d1)) ; % angle increment from previous to current location
  theta1=theta0+Dtheta ; % angle of current location
  xy(:,i)=[d1*cos(theta1)
           d1*sin(theta1)] ; 
end
xy0=xy+[Lx/2 ; Ly/2]*ones(1,M) ; % adjust to the frame center 

%% Fig. 3
Fig3=figure('Position',[100 50 600 600],'Color',[1 1 1]) ;
sft=0.07 ; 
% show emitters
Fig3a=subplot(2,2,1) ;
show8bNanoscopyImage(xy0,Lx,Ly,1,1,7,'Yes','gray','No') ; hold on
plot([100 400],Ly-[100 100],'w-',[100 100],Ly-[100-30 100+30],'w-',[400 400],Ly-[100-30 100+30],'w-') ; 
hold off
text(135,Ly-160,'300 nm','Color','white','FontSize',8)
text(100,100,'(a)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig3a,'OuterPosition',[-sft,0.5-(1-0.25)*sft+0.005,0.5+1.25*sft,0.5+1.25*sft]); 

% show far field image of all emitters 
V=CCDimage2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xy0) ;
Fig3b=subplot(2,2,2) ;
show8bimage(V,'Yes','gray','yes') ; 
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig3b,'OuterPosition',[0.5-(1-0.25)*sft+0.0115,0.5-(1-0.25)*sft+0.005,0.5+1.25*sft,0.5+1.25*sft]); 

%% Localization by a CRLB estimator
% Markov chain state space = {0, 1, 2, 3, 4} 
Nape=30 ;   % average number of activations per emitter in data movie
J=4 ;       % Maximum state
r00=0.92 ;  % Change r00 to obtain different activation Markov chains!
r01=0.5 ;   r02=0.7 ;   r03=0.8 ;   r04=1.0 ; 
r10=1-r00 ; r21=1-r01 ; r32=1-r02 ; r43=1-r03 ; 
R=[r00 r01 r02 r03 r04 % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;  % probability of de-activation
N=fix(Nape/(1-p0)) ; % total number of frames in data movie
rand('state',0) ; % initialization of pseudorandom number generator
randn('state',0) ; 
c=zeros(M,N+1) ; % states of Markov chains in data movie
for n=2:N+1, 
  for m=1:M, 
    c(m,n)=(c(m,n-1)+1)*(1-(rand<R(1,c(m,n-1)+1))) ; % state transitions
  end
end

% Estimated locations by CRLB 
a=(c~=0) ;    % a(m,n)=1 if activated; a(m,n)=0 otherwise  
Na=sum(a) ;   % number of activated emitters in nth frame
MNa=sum(Na) ; % number of emitter activations in data movie
xyActive=zeros(2,max(Na),N+1) ; % activated emitter locations in nth frame
xyCRLB=zeros(2,max(Na),N+1) ;   % locations estimated by a CRLB estimator in nth frame
xyCRLBerror=zeros(2,MNa) ;      % all estimation error vectors
xyCRLBall=zeros(2,MNa) ;        % all estimated locations 
for n=2:N+1, 
  p=0 ; 
  if Na(n)>0, 
    for m=1:M,
      if a(m,n),
        p=p+1 ;
        xyActive(:,p,n)=xy0(:,m) ;
      end
    end
  end
  if p>0,   % p=Na(n)
    [SDV SDVavg SDV_FWHM SDV_FWHMavg F_]=CRLB2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xyActive(:,1:p,n)) ;
    n1=sum(Na(1:n-1)) ; % total number of activated emitters before nth frame 
    [U L]=eig(F_) ; % eigendecomposition of inverse of Fisher information matrix
    W=U*L.^0.5*randn(2*p,1) ; % an error vector achieving Fisher information 
    for i=1:p,
      error=W(2*(i-1)+1:2*(i-1)+2) ; % error vector for ith emitter
      xyCRLB(:,i,n)=xyActive(:,i,n)+error ; % location estimated by CRLB
      xyCRLBerror(:,n1+i)=error ;
      xyCRLBall(:,n1+i)=xyCRLB(:,i,n) ;
    end
  end
  fprintf(1,'n=%3d  p=%d \n',n,p) ;
end

% Calculate RMSMD, CRLBxy, and accuracy, JAC, precision, recall
[RMSMD0 RMSMD2]=RMSMD(xy0,xyCRLBall) ;
[accuracy JAC precision recall]=Accuracy_JAC_Precesion_Recall_Core(xy0,xyCRLBall,FWHM) ;
STDx=std(xyCRLBerror(1,:)) ; % STD CRLB in x
STDy=std(xyCRLBerror(2,:)) ; % STD CRLB in y
RMSEcrlb=sqrt(STDx^2+STDy^2) ; 
fprintf(1,'r00 N 1-p0 (1-p0)*M (1-p0)*M/(Lx*Ly/10^6) (1-p0)*N \n') ; 
fprintf(1,'%5.3f %4d %5.4f %4.1f %4.2f %3.1f\n',r00,N,1-p0,(1-p0)*M,(1-p0)*M/(Lx*Ly/10^6),(1-p0)*N) ; 
fprintf(1,'r00 N RMSMD STDx STDy RMSEcrlb accuracy JAC precision recall\n') ; 
fprintf(1,'%5.3f %4d %6.1f %6.1f %6.1f %6.1f %4.1f %3.0f %3.0f %3.0f\n', ...
           r00,N,RMSMD0,STDx,STDy,RMSEcrlb,accuracy,100*JAC,100*precision,100*recall) ; 

%% 
% Show activated emitters in nth frame
% r00=0.92 considered 
Fig3b=subplot(2,2,2) ;
%n=N+1 ; 
p=sum(a(:,n)) ; 
show8bNanoscopyImage(xyActive(:,1:p,n),Lx,Ly,1,1,7,'Yes','gray','No') ; 
text(100,100,'(b)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig3b,'OuterPosition',[0.5-(1-0.25)*sft+0.0115,0.5-(1-0.25)*sft+0.005,0.5+1.25*sft,0.5+1.25*sft]); 
% show 95%-probability regions of CRLB-achieving estimators 
hold on
R095=sqrt(-2*log(1-0.95)) ; % size of ellipse for the 95%-probability region
for i=1:p,
  C=F_(2*(i-1)+1:2*(i-1)+2,2*(i-1)+1:2*(i-1)+2) ; 
  [Ui Li]=eig(C) ; 
  Wi=R095*Ui*Li.^0.5 ; % By R095: 95% estimated locations will be within the ellipse
  d=fix(400*sqrt((Li(1,1)+Li(2,2))/2)) ; % # of points on ellipse depends on variances
  theta=2*pi*(0:1/d:1) ;
  vi=[cos(theta) ; sin(theta)] ; % vi is on the unit circle
  xy_CRLBi=Wi*vi ; 
  plot(xyActive(1,i,n)+xy_CRLBi(1,:),xyActive(2,i,n)+xy_CRLBi(2,:),'r.','MarkerSize',4) ; hold on
end
hold off

% Show nth data frame 
% r00=0.92 considered 
Fig3c=subplot(2,2,3) ;
V=CCDimage2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xyActive(:,1:p,n)) ;
show8bimage(V,'Yes','gray','yes') ; 
text(15,15,'(c)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig3c,'OuterPosition',[-sft,-sft+0.01,0.5+1.25*sft,0.5+1.25*sft]); 

% Show CRLB-estimated emitters in nth frame
% r00=0.92 considered 
Fig3d=subplot(2,2,4) ;
n=N+1 ; 
p=sum(a(:,n)) ; 
show8bNanoscopyImage(xyCRLB(:,1:p,n),Lx,Ly,1,1,7,'Yes','gray','No') ; 
text(100,100,'(d)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig3d,'OuterPosition',[0.5-(1-0.25)*sft+0.0115,-sft+0.01,0.5+1.25*sft,0.5+1.25*sft]); 
%print -dtiff Fig3 % (to keep physical size, in Figure window, use File -> Save As)

%% Fig. 4 Show CRLB-estimated nanoscopy image (including all estimated emitters) 
% change r00 accordingly
Fig4=figure('Position',[100 50 600 600],'Color',[1 1 1]) ;
sft=0.07 ; 
% show Fig. 4 (a)
r00=0.98 ;  % (1-p0)*M=8.1 
r01=0.5 ;   r02=0.7 ;   r03=0.8 ;   r04=1.0 ; 
r10=1-r00 ; r21=1-r01 ; r32=1-r02 ; r43=1-r03 ;  
R=[r00 r01 r02 r03 r04 % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;  % probability of de-activation
N=fix(Nape/(1-p0)) ; % total number of frames in data movie
rand('state',0) ; % initialization of pseudorandom number generator
randn('state',0) ; 
c=zeros(M,N+1) ; % states of Markov chains in data movie
for n=2:N+1, 
  for m=1:M, 
    c(m,n)=(c(m,n-1)+1)*(1-(rand<R(1,c(m,n-1)+1))) ; % state transitions
  end
end
% Estimated locations by CRLB
a=(c~=0) ;    % a(m,n)=1 if activated; a(m,n)=0 otherwise  
Na=sum(a) ;   % number of activated emitters in nth frame
MNa=sum(Na) ; % total number of activated emitters in data movie
xyActive=zeros(2,max(Na),N+1) ;  % activated emitter locations in nth frame
xyCRLB=zeros(2,max(Na),N+1) ;    % locations estimated by a CRLB estimator in nth frame
xyCRLBerror=zeros(2,MNa) ; % all estimation error vectors
xyCRLBall=zeros(2,MNa) ;   % all estimated locations 
for n=2:N+1, 
  p=0 ; 
  if Na(n)>0, 
    for m=1:M,
      if a(m,n),
        p=p+1 ;
        xyActive(:,p,n)=xy0(:,m) ;
      end
    end
  end
  if p>0,   % p=Na(n)
    [SDV SDVavg SDV_FWHM SDV_FWHMavg F_]=CRLB2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xyActive(:,1:p,n)) ;
    n1=sum(Na(1:n-1)) ; % total number of activated emitters before nth frame 
    [U L]=eig(F_) ; % eigendecomposition of inverse of Fisher information matrix
    W=U*L.^0.5*randn(2*p,1) ; % an error vector achieving Fisher information 
    for i=1:p,
      error=W(2*(i-1)+1:2*(i-1)+2) ; % error vector for ith emitter
      xyCRLB(:,i,n)=xyActive(:,i,n)+error ; % location estimated by CRLB
      xyCRLBerror(:,n1+i)=error ;
      xyCRLBall(:,n1+i)=xyCRLB(:,i,n) ;
    end
  end
  fprintf(1,'n=%3d  p=%d \n',n,p) ;
end
%
Fig4a=subplot(2,2,1) ;
show8bNanoscopyImage(xyCRLBall,Lx,Ly,1,1,7,'Yes','gray','No') ; 
hold on
plot([100 400],Ly-[100 100],'w-',[100 100],Ly-[100-30 100+30],'w-',[400 400],Ly-[100-30 100+30],'w-') ; 
hold off
text(135,Ly-160,'300 nm','Color','white','FontSize',8)
text(100,100,'(a)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig4a,'OuterPosition',[-sft,0.5-(1-0.25)*sft+0.005,0.5+1.25*sft,0.5+1.25*sft]); 
% r00=0.98 (1-p0)*M=8.1 

% show Fig. 4 (b)
r00=0.96 ;  % (1-p0)*M=15.7 
r01=0.5 ;   r02=0.7 ;   r03=0.8 ;   r04=1.0 ; 
r10=1-r00 ; r21=1-r01 ; r32=1-r02 ; r43=1-r03 ;  
R=[r00 r01 r02 r03 r04 % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;  % probability of de-activation
N=fix(Nape/(1-p0)) ; % total number of frames in data movie
rand('state',0) ; % initialization of pseudorandom number generator
randn('state',0) ; 
c=zeros(M,N+1) ; % states of Markov chains in data movie
for n=2:N+1, 
  for m=1:M, 
    c(m,n)=(c(m,n-1)+1)*(1-(rand<R(1,c(m,n-1)+1))) ; % state transitions
  end
end
% Estimated locations by CRLB
a=(c~=0) ;    % a(m,n)=1 if activated; a(m,n)=0 otherwise  
Na=sum(a) ;   % number of activated emitters in nth frame
MNa=sum(Na) ; % total number of activated emitters in data movie
xyActive=zeros(2,max(Na),N+1) ;  % activated emitter locations in nth frame
xyCRLB=zeros(2,max(Na),N+1) ;    % locations estimated by a CRLB estimator in nth frame
xyCRLBerror=zeros(2,MNa) ; % all estimation error vectors
xyCRLBall=zeros(2,MNa) ;   % all estimated locations 
for n=2:N+1, 
  p=0 ; 
  if Na(n)>0, 
    for m=1:M,
      if a(m,n),
        p=p+1 ;
        xyActive(:,p,n)=xy0(:,m) ;
      end
    end
  end
  if p>0,   % p=Na(n)
    [SDV SDVavg SDV_FWHM SDV_FWHMavg F_]=CRLB2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xyActive(:,1:p,n)) ;
    n1=sum(Na(1:n-1)) ; % total number of activated emitters before nth frame 
    [U L]=eig(F_) ; % eigendecomposition of inverse of Fisher information matrix
    W=U*L.^0.5*randn(2*p,1) ; % an error vector achieving Fisher information 
    for i=1:p,
      error=W(2*(i-1)+1:2*(i-1)+2) ; % error vector for ith emitter
      xyCRLB(:,i,n)=xyActive(:,i,n)+error ; % location estimated by CRLB
      xyCRLBerror(:,n1+i)=error ;
      xyCRLBall(:,n1+i)=xyCRLB(:,i,n) ;
    end
  end
  fprintf(1,'n=%3d  p=%d \n',n,p) ;
end
%
Fig4b=subplot(2,2,2) ;
show8bNanoscopyImage(xyCRLBall,Lx,Ly,1,1,7,'Yes','gray','No') ; 
hold off
text(100,100,'(b)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig4b,'OuterPosition',[0.5-(1-0.25)*sft+0.0115,0.5-(1-0.25)*sft+0.005,0.5+1.25*sft,0.5+1.25*sft]); 
% r00=0.96 (1-p0)*M=15.7

% show Fig. 4 (c)
r00=0.94 ;  % (1-p0)*M=22.9 
r01=0.5 ;   r02=0.7 ;   r03=0.8 ;   r04=1.0 ; 
r10=1-r00 ; r21=1-r01 ; r32=1-r02 ; r43=1-r03 ;  
R=[r00 r01 r02 r03 r04 % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;  % probability of de-activation
N=fix(Nape/(1-p0)) ; % total number of frames in data movie
rand('state',0) ; % initialization of pseudorandom number generator
randn('state',0) ; 
c=zeros(M,N+1) ; % states of Markov chains in data movie
for n=2:N+1, 
  for m=1:M, 
    c(m,n)=(c(m,n-1)+1)*(1-(rand<R(1,c(m,n-1)+1))) ; % state transitions
  end
end
% Estimated locations by CRLB
a=(c~=0) ;    % a(m,n)=1 if activated; a(m,n)=0 otherwise  
Na=sum(a) ;   % number of activated emitters in nth frame
MNa=sum(Na) ; % total number of activated emitters in data movie
xyActive=zeros(2,max(Na),N+1) ;  % activated emitter locations in nth frame
xyCRLB=zeros(2,max(Na),N+1) ;    % locations estimated by a CRLB estimator in nth frame
xyCRLBerror=zeros(2,MNa) ; % all estimation error vectors
xyCRLBall=zeros(2,MNa) ;   % all estimated locations 
for n=2:N+1, 
  p=0 ; 
  if Na(n)>0, 
    for m=1:M,
      if a(m,n),
        p=p+1 ;
        xyActive(:,p,n)=xy0(:,m) ;
      end
    end
  end
  if p>0,   % p=Na(n)
    [SDV SDVavg SDV_FWHM SDV_FWHMavg F_]=CRLB2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xyActive(:,1:p,n)) ;
    n1=sum(Na(1:n-1)) ; % total number of activated emitters before nth frame 
    [U L]=eig(F_) ; % eigendecomposition of inverse of Fisher information matrix
    W=U*L.^0.5*randn(2*p,1) ; % an error vector achieving Fisher information 
    for i=1:p,
      error=W(2*(i-1)+1:2*(i-1)+2) ; % error vector for ith emitter
      xyCRLB(:,i,n)=xyActive(:,i,n)+error ; % location estimated by CRLB
      xyCRLBerror(:,n1+i)=error ;
      xyCRLBall(:,n1+i)=xyCRLB(:,i,n) ;
    end
  end
  fprintf(1,'n=%3d  p=%d \n',n,p) ;
end
%
Fig4c=subplot(2,2,3) ;
show8bNanoscopyImage(xyCRLBall,Lx,Ly,1,1,7,'Yes','gray','No') ; 
hold off
text(100,100,'(c)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig4c,'OuterPosition',[-sft,-sft+0.01,0.5+1.25*sft,0.5+1.25*sft]); 
% r00=0.94 (1-p0)*M=22.9 

% show Fig. 4 (d)
r00=0.92 ;  % (1-p0)*M=29.6
r01=0.5 ;   r02=0.7 ;   r03=0.8 ;   r04=1.0 ; 
r10=1-r00 ; r21=1-r01 ; r32=1-r02 ; r43=1-r03 ;  
R=[r00 r01 r02 r03 r04 % matrix of state transition probabilities
   r10 0   0   0   0
   0   r21 0   0   0
   0   0   r32 0   0
   0   0   0   r43 0] ;
den=1+r10+r10*r21+r10*r21*r32+r10*r21*r32*r43 ; 
p0=1/den ;  % probability of de-activation
N=fix(Nape/(1-p0)) ; % total number of frames in data movie
rand('state',0) ; % initialization of pseudorandom number generator
randn('state',0) ; 
c=zeros(M,N+1) ; % states of Markov chains in data movie
for n=2:N+1, 
  for m=1:M, 
    c(m,n)=(c(m,n-1)+1)*(1-(rand<R(1,c(m,n-1)+1))) ; % state transitions
  end
end
% Estimated locations by CRLB
a=(c~=0) ;    % a(m,n)=1 if activated; a(m,n)=0 otherwise  
Na=sum(a) ;   % number of activated emitters in nth frame
MNa=sum(Na) ; % total number of activated emitters in data movie
xyActive=zeros(2,max(Na),N+1) ;  % activated emitter locations in nth frame
xyCRLB=zeros(2,max(Na),N+1) ;    % locations estimated by a CRLB estimator in nth frame
xyCRLBerror=zeros(2,MNa) ; % all estimation error vectors
xyCRLBall=zeros(2,MNa) ;   % all estimated locations 
for n=2:N+1, 
  p=0 ; 
  if Na(n)>0, 
    for m=1:M,
      if a(m,n),
        p=p+1 ;
        xyActive(:,p,n)=xy0(:,m) ;
      end
    end
  end
  if p>0,   % p=Na(n)
    [SDV SDVavg SDV_FWHM SDV_FWHMavg F_]=CRLB2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xyActive(:,1:p,n)) ;
    n1=sum(Na(1:n-1)) ; % total number of activated emitters before nth frame 
    [U L]=eig(F_) ; % eigendecomposition of inverse of Fisher information matrix
    W=U*L.^0.5*randn(2*p,1) ; % an error vector achieving Fisher information 
    for i=1:p,
      error=W(2*(i-1)+1:2*(i-1)+2) ; % error vector for ith emitter
      xyCRLB(:,i,n)=xyActive(:,i,n)+error ; % location estimated by CRLB
      xyCRLBerror(:,n1+i)=error ;
      xyCRLBall(:,n1+i)=xyCRLB(:,i,n) ;
    end
  end
  fprintf(1,'n=%3d  p=%d \n',n,p) ;
end
%
Fig4d=subplot(2,2,4) ;
show8bNanoscopyImage(xyCRLBall,Lx,Ly,1,1,7,'Yes','gray','No') ; 
hold off
text(100,100,'(d)','Color','white')
set(gca,'XTick',[]) ; % Turn off X and Y ticks 
set(gca,'YTick',[]) ; 
axis off
set(Fig4d,'OuterPosition',[0.5-(1-0.25)*sft+0.0115,-sft+0.01,0.5+1.25*sft,0.5+1.25*sft]); 
% r00=0.92 (1-p0)*M=29.6
% print -dtiff Fig4 % (to keep physical size, in Figure window, use File -> Save As)

%% Fig. 5 Show metrics 
% Result and figures
% Nape=30 ;Ih=300000 ; SPNR=0.2 ; SGNR=0.3 ; 
%    r00   N  1-p0 (1-p0)*M (1-p0)*M/(Lx*Ly) (1-p0)*N 
Np=[0.995 3601 0.0083  2.1 0.50 30.0
    0.99  1815 0.0165  4.1 0.98 30.0
    0.98   922 0.0325  8.1 1.94 30.0
    0.97   625 0.0480 12.0 2.86 30.0
    0.96   476 0.0630 15.7 3.75 30.0
    0.95   387 0.0775 19.4 4.62 30.0
    0.94   327 0.0916 22.9 5.46 29.9
    0.93   285 0.1052 26.3 6.27 30.0
    0.92   253 0.1185 29.6 7.06 30.0
    0.91   228 0.1313 32.8 7.83 29.9] ;
%    r00    N    RMSMD    STDx   STDy  RMSEcrlb accuracy JAC precision recall
Mt=[0.995 3601     4.7     3.5    3.5     5.0  4.8 100 100 100
    0.990 1815     6.0     4.5    4.6     6.4  6.1 100 100 100
    0.980  922    11.2    11.1    9.8    14.8 10.0 100 100 100
    0.970  625    22.9    25.5   26.6    36.8 14.1 100 100 100
    0.960  476    38.9    38.1   50.1    63.0 17.7 100 100 100
    0.950  387    57.7    74.9   58.4    95.0 22.7  99  99 100
    0.940  327   132.2   106.1  150.1   183.8 24.9  99  99 100
    0.930  285   705.9   646.3  423.9   772.9 28.9  97  97 100
    0.920  253  2515.9  2129.1 1503.7  2606.6 35.6  96  96 100
    0.910  228 13636.9 10214.6 9407.8 13886.9 39.2  93  93 100] ;
  
Fig5=figure('Position',[100 100 620 300],'Color',[1 1 1]) ;

% Fig5. (a) RMSMD, and accuracy 
Fig5a=subplot(1,2,1) ;
clear lg ; 
lg(1)=plot(Np(:,4),Mt(:,3),'k-*','LineWidth',0.125,'MarkerSize',5) ; hold on
lg(2)=plot(Np(:,4),Mt(:,7),'k-o','LineWidth',0.125,'MarkerSize',5) ; 
abcd=[8.1 15.7 22.9 29.6] ; 
for i=1:4,
  plot(abcd(i)*[1 1],[0,3000],'k:','LineWidth',0.125) ; 
end
hold off
axis([5 30 0 3000]) ; 
legend(lg,'RMSMD','Accuracy','location','northwest') ;
xlabel('Ma','FontSize',8) ;
ylabel('RMSMD, Accuracy (nm)','FontSize',8) ; 
text(6,180,'(a)','Color','black')
set(gca,'FontSize',8) ; 
set(Fig5a,'OuterPosition',[0,0,0.488,1]); 

% Fig5. (b) precision, recall, JAC
Fig5b=subplot(1,2,2) ;
% precision 
lg(1)=plot(Np(:,4),Mt(:,9),'k-o','LineWidth',0.125,'MarkerSize',5) ; hold on
% recall 
lg(2)=plot(Np(:,4),Mt(:,10),'k-s','LineWidth',0.125,'MarkerSize',5) ; 
% JAC
lg(3)=plot(Np(:,4),Mt(:,9),'k-x','LineWidth',0.125,'MarkerSize',5) ; 
abcd=[8.1 15.7 22.9 29.6] ; 
for i=1:4,
  plot(abcd(i)*[1 1],[0,3000],'k:','LineWidth',0.125) ; 
end
hold off
axis([5 30 0 120]) ; 
legend(lg,'Precision','Recall','JAC','Location','southeast') ;
xlabel('Ma','FontSize',8) ;
ylabel('Precision, Recall, JAC (%)','FontSize',8) ; 
text(6,7,'(b)','Color','black')
set(gca,'FontSize',8) ; 
set(Fig5b,'OuterPosition',[0.51,0,0.515,1]); 
% print -desp Fig5 (to keep physical size, in Figure window, use File -> Save As)

%% Show sequence of estimated nanoscopy images 
figure('Position',[100 50 600 600],'Color',[1 1 1]) ;
sft=0.07 ; 
FPS=10 ; % frames/s in show
for n=2:N+1,
  p=sum(a(:,n)) ; % number of activated emitters in nth frame
  fprintf(1,'N=%d n=%3d  p=%d \n',N,n,p) ;
  t0=cputime ;
%
  Figa=subplot(2,2,1) ;
  show8bNanoscopyImage(xyActive(:,1:p,n),Lx,Ly,1,1,7,'Yes','gray','No') ;
  hold on
  plot([100 400],Ly-[100 100],'w-',[100 100],Ly-[100-30 100+30],'w-',[400 400],Ly-[100-30 100+30],'w-') ;
  hold off
  text(135,Ly-160,'300 nm','Color','white','FontSize',8)
  set(gca,'XTick',[]) ; % Turn off X and Y ticks
  set(gca,'YTick',[]) ;
  axis off
  set(Figa,'OuterPosition',[-sft,0.5-(1-0.25)*sft+0.005,0.5+1.25*sft,0.5+1.25*sft]);
%
  Figb=subplot(2,2,2) ;
  if p>0,
    V=CCDimage2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xyActive(:,1:p,n)) ;
  else
    V=zeros(Ky,Kx) ;
  end
  show8bimage(V,'Yes','gray','Yes') ;
  axis off
  set(gca,'XTick',[]) ; % Turn off X and Y ticks
  set(gca,'YTick',[]) ;
  set(Figb,'OuterPosition',[0.5-(1-0.25)*sft+0.0115,0.5-(1-0.25)*sft+0.005,0.5+1.25*sft,0.5+1.25*sft]);
%
  Figc=subplot(2,2,3) ;  % show current estimated locations
  if p>0,
    show8bNanoscopyImage(xyCRLB(:,1:p,n),Lx,Ly,1,1,7,'Yes','gray','No') ;
  else
    show8bNanoscopyImage([],Lx,Ly,1,1,7,'Yes','gray','No') ;
  end
  set(gca,'XTick',[]) ; % Turn off X and Y ticks
  set(gca,'YTick',[]) ;
  axis off
  set(Figc,'OuterPosition',[-sft,-sft+0.01,0.5+1.25*sft,0.5+1.25*sft]);
%
  Figd=subplot(2,2,4) ;  % show all estimated locations
  pt=sum(sum(a(:,1:n))) ; % number of activated emitters from 1st to nth frames
  show8bNanoscopyImage(xyCRLBall(:,1:pt),Lx,Ly,1,1,7,'Yes','gray','No') ; 
  set(gca,'XTick',[]) ; % Turn off X and Y ticks
  set(gca,'YTick',[]) ;
  axis off
  set(Figd,'OuterPosition',[0.5-(1-0.25)*sft+0.0115,-sft+0.01,0.5+1.25*sft,0.5+1.25*sft]);
  gf=getframe(gcf) ;	% show image
  while cputime-t0<1/FPS, end
end

