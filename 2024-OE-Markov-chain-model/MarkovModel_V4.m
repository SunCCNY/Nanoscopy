%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates figures and data movies for the paper titled 
% "Markov chain models of emitter activations in single molecule 
% localization microscopy". 
% 
% Section(A): Cycled illumination for 3D imaging 
% Section(B): Continuously illumination for 2D imaging 
% 
% You can choose to run one section
% 
% Yi Sun
% Electrical Engineering Department
% The City College of City University of New York
% E-mail: ysun@ccny.cuny.edu
% 01/06, 04/29, 05/07, 08/15/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section (A): Cycled illumination for 3D imaging 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Emitter distance 
clear
fprintf(1,'Section (A): Cycled illumination for 3D imaging \n') ; 
eD=40 ;         % nm 
rng('default') ; 
rng(eD) ;       % use eD as key to initialize random number generators
fprintf(1,'Emitter lateral distance: %d (nm) \n',eD) ; 
%%% Optical system: 3D astigmatic PSF
c=205 ;
d=290 ;         % depth of microscope
sigmax0=280/2 ; Ax=0.05 ; Bx=0.03 ;   % nm
sigmay0=270/2 ; Ay=-0.01 ; By=0.02 ;	% 
%%% show sigmax, sigmay vs. z for 3d PSF
z=-600:600 ;
sigmax=sigmax0*(1+(z+c).^2/d^2+Ax*(z+c).^3/d^3+Bx*(z+c).^4/d^4).^0.5 ;
sigmay=sigmay0*(1+(z-c).^2/d^2+Ay*(z-c).^3/d^3+By*(z-c).^4/d^4).^0.5 ;
figure('Units','inches','Position',[4 6 4*[1 0.75]],'Color',[1 1 1]) ;
lg=plot(z,sigmax,'r-',z,sigmay,'b--','LineWidth',1) ;
axis([-600 600 0 500])
grid on
legend(lg,'\sigma_x(z)','\sigma_y(z)','Location','southeast')
xlabel('z (nm)')
ylabel('\sigma_x(z), \sigma_y(z) (nm)')
set(gca,'FontSize',8,'FontWeight','bold') ;

%%% Frame parameters
% sample is located at [0,Lx]x[0,Ly]x[-Lz,Lz]
Lx=4e3 ; Ly=Lx ;      % frame size in nm
Dx=1e2 ; Dy=1e2 ;     % pixel size in nm
Kx=Lx/Dx ; Ky=Ly/Dy ; % frame size in pixels
Lz=400 ;              % 2 x Lz = axial depth in nm
%%% Emitter intensity and signal to noise ratio
Dt=0.01 ;             % second, time per frame (1/Dt is frame rate) 
Ih=240000 ;           % average number of detected photons per emitter per second
DtIh=Dt*Ih ;          % photon count per frame per emitter 
% 'medium SNR'        % 
b=0.3 ;               % mean of Poisson noise (photons/s/nm^2)
                      % uniform background Poisson noise is assumed
G=0.2 ;               % variance of Gaussian noise (photons/s/nm^2) 
rp=Ih/b ;             % 800000 (nm^2/emitter) 
rg=Ih/G ;             % 1200000 (nm^2/emitter) 
r=rp*rg/(rp+rg) ;     % 480000 (nm^2/emitter) 
mu=5 ;                % mean of Gaussian noise (photons/s/nm^2)
Coff=mu*Dt*Dx*Dy ;    % Coff=819.2 photons/pixel; Camera offset in effect
% calculate SNR by Ref. [17] 
z0=0 ;                % consider z0=0 only 
sigmax_z0= ...        % 172.8791 nm
  sigmax0*(1+(z0+c).^2/d^2+Ax*(z0+c).^3/d^3+Bx*(z0+c).^4/d^4).^0.5 ;
sigmay_z0= ...        % 165.7935 nm
  sigmay0*(1+(z0-c).^2/d^2+Ay*(z0-c).^3/d^3+By*(z0-c).^4/d^4).^0.5 ;
SNR=10*log10(r) ...   % 1.2193 (dB) 
  -10*log10(sigmax_z0*sigmay_z0)-11.02

%%% Emitter locations - ground truth
fprintf(1,'Emitter locations: \n') ;
M=500 ;               % number of emitters
phi=0.1 ; beta=7.5+0.5*rand ;  % helix parameters 
xy1=zeros(3,M) ; 
J=zeros(1,M) ; 
m=1 ; J(m)=2 ;        % 1st emitter locaiton
xy1(:,m)=[beta*J(m)*cos(phi*J(m))  
          beta*J(m)*sin(phi*J(m))
          -Lz+(m-0.5)*(2*Lz/M)] ; 
for m=2:M
  if mod(m,10)==0||m==1
    fprintf(1,'M=%3d m=%3d \n',M,m) ;
  end
  syms t0
  St=vpasolve((beta*t0*cos(phi*t0)-xy1(1,m-1))^2+(beta*t0*sin(phi*t0)-xy1(2,m-1))^2==eD^2,t0,J(m-1)+0.1) ;
  J(m)=St ;
  xy1(:,m)=[beta*J(m)*cos(phi*J(m)) 
            beta*J(m)*sin(phi*J(m))
            -Lz+(m-0.5)*(2*Lz/M)] ; 
end
phi_=2*pi*rand ;            % random initial position in lateral plane
xy_=[cos(phi_) -sin(phi_) 
     sin(phi_) cos(phi_)]*xy1(1:2,:)+2*randn(2,1) ; 
xy0=[(max(xy_(1,:))+min(xy_(1,:)))/2 ; (max(xy_(2,:))+min(xy_(2,:)))/2] ;
xy_=xy_-xy0+[Lx/2 ; Ly/2] ; % adjust to frame center in lateral plane
xyz=[xy_ ; xy1(3,:)] ;      % ground truth emitter locaitons 

%%% Fig. 8. 3D emitter locations 
fprintf(1,'Fig. 8 \n') ; 
wx=3 ; wy=0.6*wx ; dx=0.62 ; dy=-0.12 ; 
PS=[dx/2+0.04 dy/2+0.25 wx-dx wy-dy] ;    % Positions of subfigures 
figure('Units','inches','Position',[3 4 wx wy],'Color',[1 1 1]) ;
subplot(1,1,1,'Units','inches','Position',PS) ;
plot3(xyz(1,:),xyz(2,:),xyz(3,:),'k.','MarkerSize',3) 
axis([0 Lx 0 Ly -Lz Lz])
xlabel('x (nm)','FontSize',6)
ylabel('y (nm)','FontSize',6)
zlabel('z (nm)','FontSize',6)
set(gca,'XTick',[0 1000 2000 3000 4000]) ;
set(gca,'XTickLabel',[0 1000 2000 3000 4000]) ;
set(gca,'YTick',[0 1000 2000 3000 4000]) ;
set(gca,'YTickLabel',[0 1000 2000 3000 4000]) ;
set(gca,'ZTick',[-400 -200 0 200 400]) ;
set(gca,'ZTickLabel',[-400 -200 0 200 400]) ;
set(gca,'FontSize',6) ;
print('Fig8.tif','-dtiffn')

%%% Emitter activation parameters: Cycled illumination, state 0, 1, 2 
% A-phase
tA0=0.014 ;       % < tA1, tD0, 1.44*KA*Dt=0.0144*KA, mean of TA0 in sec
pA0=exp(-Dt/tA0)  % 0.4895, probability to retain state 0 in Dt
tA1=0.02 ;        % > tA0, tD1, mean of time TA1 in sec
pA1=exp(-Dt/tA1)  % 0.6065, probability to retain state 1 in Dt
% D-phase 
tD0=0.30 ;        % > tA0, tD1, mean of time TD0 in sec 
pD0=exp(-Dt/tD0)  % 0.9672, probability to retain state 0 in Dt
tD1=0.008 ;       % < tA1, tD0, 1.44*KD*Dt=0.0144*KD, mean of time TD1 in sec
pD1=exp(-Dt/tD1)  % 0.2865, probability to retain state 1 in Dt
% photoactivatable period 
t=5 ;             % mean of time T 
p=exp(-Dt/t)      % 0.9980, probability to retain photoactivatable in Dt
% transition probabilities
a00=pA0*p         % 0.4886, Trans. prob. from state 0 to 0
a10=(1-pA0)*p     % 0.5094, Trans. prob. from state 0 to 1
a20=1-p           % 0.0020, Trans. prob. from state 0 to 2
a01=(1-pA1)*p     % 0.3927, Trans. prob. from state 1 to 0
a11=pA1*p         % 0.6053, Trans. prob. from state 1 to 1
a21=1-p           % 0.0020, Trans. prob. from state 1 to 2
d00=pD0*p         % 0.9653, Trans. prob. from state 0 to 0
d10=(1-pD0)*p     % 0.0327, Trans. prob. from state 0 to 1
d20=1-p           % 0.0020, Trans. prob. from state 0 to 2
d01=(1-pD1)*p     % 0.7121, Trans. prob. from state 1 to 0
d11=pD1*p         % 0.2859, Trans. prob. from state 1 to 1
d21=1-p           % 0.0020, Trans. prob. from state 1 to 2

%%% Emitter activation states: KA=1, KD=3 
C=125 ;                     % number of cycles
KA=1 ;                      % number of A-frames per cycle
KD=3 ;                      % number of D-frames per cycle
K=KA+KD ;                   % total number of frames per cycle 
N=C*K ;                     % number of frames in data movie
                            % Temporal resolution (TR): N*Dt=5 sec
[ca,p,h1]=emActMarkovCycle(t,tA1,tA0,tD1,tD0,Dt,C,KA,KD,M) ; 
a=(ca==1) ;                 % a(m,n)=1 if activated; a(m,n)=0 otherwise  
                            % sum(sum(a')==0): # of emitters never activated 
Ma_=sum(a) ;                % number of activated emitters in nth frame
% compare formulas and estimates 
Np=p*(1-p^N)/(1-p)          % 315.74
Np_=sum(sum(ca~=2))/M       % 313.64
Nae1=sum(p.^(1:K).*h1)      % 0.80
Nae=((1-p^N)/(1-p^K))*Nae1  % 63.65
Nae_=sum(Ma_)/M             % 63.04
Na=M*Nae                    % 31826
Na_=M*Nae_                  % 31520

%%% Fig. 9. Number of activated emitters Ma(n) VS frame n
fprintf(1,'Fig. 9 \n') ; 
str=['r-' ; 'g-' ; 'm-' ; 'b-' ; 'c-'] ; 
wx=2.6 ; wy=0.75*wx ; dx=0.4 ; dy=0.34 ; 
PS=[dx/2+0.1 dy/2+0.1 wx-dx wy-dy] ;    % Positions of subplot 
figure('Units','inches','Position',[3 4 wx wy],'Color',[1 1 1]) ;
subplot(1,1,1,'Units','inches','Position',PS) ;
plot((1:N),Ma_,'k.','MarkerSize',5) ; hold on 
lg=zeros(1,K) ; 
for j=1:K
  n=(0:C-1)*K+j ; 
  Ma=M*p.^n*h1(j) ;        % average # A-frames for M emitters in nth frame
  lg(j)=plot(n,Ma,str(j,:),'MarkerSize',5) ; 
end
hold off 
legend(lg,'Frame 1','Frame 2','Frame 3','Frame 4', ...
  'Location','northeast','FontSize',6) ; 
axis([1 N 0 300])
ylabel('Number of activated emitters','FontSize',6)
xlabel('Frame n','FontSize',6) 
set(gca,'XTick',[0 50 100 150 200 250 300 350 400 450 500]) ;
set(gca,'XTickLabel',[0 50 100 150 200 250 300 350 400 450 500]) ;
set(gca,'YTick',[0 50 100 150 200 250 300]) ;
set(gca,'YTickLabel',[0 50 100 150 200 250 300]) ;
set(gca,'FontSize',6) ;
grid on
print('Fig9.tif','-dtiffn')

%%% Stack activated emitters in each frame
xyzActive=zeros(3,M,N) ;  % activated emitter locations in nth frame
for n=1:N
  if Ma_(n)>0
    k=0 ; 
    for m=1:M
      if a(m,n)
        k=k+1 ;
        xyzActive(:,k,n)=xyz(:,m) ;
      end
      if k==Ma_(n), break ; end
    end
  end
end

%%% Fig. 10: Data frames 101, 102, 103, 104
fprintf(1,'Fig. 10 \n') ; 
n12=101:104 ; 
em=[11 105 205 339 409] ;       % indices of emitters to be circled
st=['go' ; 'mo' ; 'yo' ; 'co' ; 'bo'] ;
wx=1.3 ; wy=1.3 ; dx=0.04 ; dy=0.04 ; 
figure('Units','inches','Position',[1 2 4*wx wy],'Color',[1 1 1]) ;
PS=[0*wx+dx/2 dy/2 wx-dx wy-dy  % Positions of subplots 
    1*wx+dx/2 dy/2 wx-dx wy-dy
    2*wx+dx/2 dy/2 wx-dx wy-dy
    3*wx+dx/2 dy/2 wx-dx wy-dy] ;
for i=1:length(n12)
  n=n12(i) ; 
  fprintf(1,'n=%d  Ma_(n)=%d \n',n,Ma_(n)) ; 
  k=Ma_(n) ;
  if k>0
    U=AS3D_Frame(c,d,sigmax0,Ax,Bx,sigmay0,Ay,By,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,xyzActive(:,1:k,n)) ;
  else
    U=zeros(Ky,Kx) ;
  end
  subplot(1,4,i,'Units','inches','Position',PS(i,:)) ;
  show8bimage(U,'Yes','gray','No') ; hold on
  plot(xyz(1,:)/Dx+0.5,xyz(2,:)/Dy+0.5,'w.',xyzActive(1,1:Ma_(n),n)/Dx+0.5, ...
    xyzActive(2,1:Ma_(n),n)/Dy+0.5,'r.','MarkerSize',3.5) ;
  for m=1:length(em)
    plot(xyz(1,em(m))/Dx+0.5,xyz(2,em(m))/Dy+0.5,st(m,:),'MarkerSize',3.5) ;
  end
  if i==1
      plot([2 12],Ky-[1 1],'w-',[2 2],Ky-[1-0.8 1+0.8],'w-',[12 12],Ky-[1-0.8 1+0.8],'w-') ;
      text(3.5,Ky-3.5,'1 {\mu}m','Color','white','FontSize',6) ;
  end
  hold off
  fm=compose('Frame %d',n) ; 
  text(27,3,fm,'Color','white','fontsize',6) ;
  axis off
  [min(min(U)) max(max(U))]
end
print('Fig10.tif','-dtiffn')
[[n12 ; a(:,n12)] [0 ; sum(a(:,n12),2)] (0:N)']   
Ma_(n12)   % [229    74    30    21]

%%% Generate a video of data movie: 100 frames  
fprintf(1,'Generate a video of data movie: 100 frames\n') ; 
vidObj=VideoWriter('Cycled3D.avi') ; 
vidObj.FrameRate=5 ;      % frames/s in movie
open(vidObj) ; 
wx=3.4 ; wy=wx ; dx=0.1 ; dy=0.1 ; 
PS=[dx/2 dy/2 wx-dx wy-dy] ;    % Positions of subfigures 
figure('Units','inches','Position',[3 4 wx wy],'Color',[1 1 1]) ;
subplot(1,1,1,'Units','inches','Position',PS) ;
for n=1:100
  k=Ma_(n) ;
  if k>0
    U=AS3D_Frame(c,d,sigmax0,Ax,Bx,sigmay0,Ay,By,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,xyzActive(:,1:k,n)) ;
  else
    U=zeros(Ky,Kx) ;
  end
  show8bimage(U,'Yes','gray','No') ; hold on
  plot(xyz(1,:)/Dx+0.5,xyz(2,:)/Dy+0.5,'w.',xyzActive(1,1:k,n)/Dx+0.5,xyzActive(2,1:k,n)/Dy+0.5,'r.') ;
  plot([2 7],Ky-[1 1],'w-',[2 2],Ky-[1-0.5 1+0.5],'w-',[7 7],Ky-[1-0.5 1+0.5],'w-') ;
  hold off
  text(2.2,Ky-2.3,'500 nm','Color','white','FontSize',8) ;
  currFrame=getframe(gcf) ;
  writeVideo(vidObj,currFrame) ;
  fprintf(1,'N=%d n=%3d Ma_(n)=%3d \n',N,n,k) ;
end
close(vidObj) ;

%%% Generate and save a data movie
fprintf(1,'Generate a data movie: \n') ; 
for n=1:N
  k=Ma_(n) ;
  if k>0
    U=AS3D_Frame(c,d,sigmax0,Ax,Bx,sigmay0,Ay,By,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,xyzActive(:,1:k,n)) ; 
  else
    U=zeros(Ky,Kx) ; 
  end
  U16=uint16(U) ;
  filename_Frame=strcat('Movie-Cycled3D\Cycled3D-',num2str(n),'.tif') ;
  imwrite(U16,filename_Frame) ;    % save data frame
  if mod(n,10)==0||n==1
    fprintf(1,'N=%3d n=%3d Ma_(n)=%d \n',N,n,Ma_(n)) ;
  end
end

%%% Verify saved frames: 101, 102, 103, 104
fprintf(1,'Verification of frames 101-104 \n') ; 
n12=101:104 ; 
em=[11 105 205 339 409] ;           % show emitters in em
st=['go' ; 'mo' ; 'yo' ; 'co' ; 'bo'] ; 
wx=1.7 ; wy=1.7 ; dx=0.04 ; dy=0.04 ; 
figure('Units','inches','Position',[1 2 4*wx wy],'Color',[1 1 1]) ;
PS=[0*wx+dx/2 dy/2 wx-dx wy-dy      % Positions of subfigures 
    1*wx+dx/2 dy/2 wx-dx wy-dy
    2*wx+dx/2 dy/2 wx-dx wy-dy
    3*wx+dx/2 dy/2 wx-dx wy-dy] ;
for i=1:4 %length(n12)
  n=n12(i) ; 
  fprintf(1,'n=%d  Ma_(n)=%d \n',n,Ma_(n)) ; 
  filename_Frame=strcat('Movie-Cycled3D\Cycled3D-',num2str(n),'.tif') ;
  U16=imread(filename_Frame); % read a frame
  U=double(U16) ;
  subplot(1,4,i,'Units','inches','Position',PS(i,:)) ;
  show8bimage(U,'Yes','gray','No') ; hold on
  plot(xyz(1,:)/Dx+0.5,xyz(2,:)/Dy+0.5,'w.',xyzActive(1,1:Ma_(n),n)/Dx+0.5, ...
    xyzActive(2,1:Ma_(n),n)/Dy+0.5,'r.','MarkerSize',3.5) ;
  for m=1:length(em)
    plot(xyz(1,em(m))/Dx+0.5,xyz(2,em(m))/Dy+0.5,st(m,:),'MarkerSize',3.5) ;
  end
  if i==1
      plot([2 7],Ky-[1 1],'w-',[2 2],Ky-[1-0.5 1+0.5],'w-',[7 7],Ky-[1-0.5 1+0.5],'w-') ;
      text(1.8,Ky-2.5,'500 nm','Color','white','FontSize',5) ;
  end
  hold off
  fm=compose('Frame %d',n) ; 
  text(29.5,3,fm,'Color','white','fontsize',6) ;
  axis off
  [min(min(U)) max(max(U))]
end
[[n12 ; a(:,n12)] [0 ; sum(a(:,n12),2)] (0:N)']   
Ma_(n12)   % [229    74    30    21]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section (B): Continuous illumination for 2D imaging 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Emitter distance 
clear
fprintf(1,'Section (B): Continuous illumination for 2D imaging \n') ; 
eD=40 ;               % nm 
rng('default') ; 
rng(eD) ;             % use eD as key to initialize random number generators
fprintf(1,'Emitter lateral distance: %d (nm) \n',eD) ; 
%%% Optical system 
na=1.40 ; 
lambda=723 ;                  % Alexa700 wavelength in nm
a=2*pi*na/lambda ;  
% 2D Gaussian PSF; sigma is estimated from Airy PSF
sigma=1.3238/a ;              % sigma=108.81; 2*sigma=217.61 (nm) 
FWHM=2*sqrt(2*log(2))*sigma ; % FWHM=256.22 (nm)
%%% Frame paramiters
% Region of view: [0,Lx]x[0,Ly]
Lx=4e3 ; Ly=Lx ;      % frame size in nm
Dx=1e2 ; Dy=1e2 ;     % pixel size of cammera
Kx=Lx/Dx ; Ky=Ly/Dy ; % frame size in pixels
%%% Emitter intensity and signal to noise ratio
Dt=0.01 ;             % second, time per frame (1/Dt is frame rate) 
Ih=300000 ;           % average number of detected photons per emitter per second
DtIh=Dt*Ih ;          % photon count per frame per emitter 
% 'mediumSNR'         % 
b=5 ;                 % mean of Poisson noise (photons/s/nm^2)
                      % uniform backgroun Poisson noise is assumed
G=3 ;                 % variance of Gaussian noise (photons/s/nm^2)
% calcaulate SNR by Ref. [18]
betas=0.07912 ;       % Ref. [18]
beta=betas/sigma^2 ; 
rp=Ih/b ;             % 60000
nup=beta*rp ; 
SPNR=10*log10(nup)    % -3.97 (dB)
rg=Ih/G ;             % 100000
nug=beta*rg ; 
SGNR=10*log10(nug)    % -1.75 (dB)
r=rp*rg/(rp+rg) ;     % 37500
nu=beta*r ; 
SNR=10*log10(nu)      % -6.01 (dB): SNR=10*log10(r)-20*log10(sigma)-11.02
mu=5 ;                % mean of Gaussian noise (photons/s/nm^2)
Coff=mu*Dt*Dx*Dy ;    % Coff=819.2 photons/pixel; Camera offset in effect

%%% Emitter locations - ground truth
fprintf(1,'Emitter locations: \n') ;
M=500 ;
phi=0.1 ; beta=7+0.5*rand ;  % helix parameters xy1=zeros(2,M) ;
xy1=zeros(2,M) ; 
J=zeros(1,M) ; 
m=1 ; J(m)=2 ;            % 1st emitter locaiton
xy1(:,m)=[beta*J(m)*cos(phi*J(m)) ; beta*J(m)*sin(phi*J(m))] ; 
for m=2:M
  if mod(m,50)==0||m==2
    fprintf(1,'M=%3d m=%3d \n',M,m) ;
  end
  syms t0
  St=vpasolve((beta*t0*cos(phi*t0)-xy1(1,m-1))^2+(beta*t0*sin(phi*t0)-xy1(2,m-1))^2==eD^2,t0,J(m-1)+0.1) ;
  J(m)=St ;
  xy1(:,m)=[beta*J(m)*cos(phi*J(m)) ; beta*J(m)*sin(phi*J(m))] ;
end
phi_=2*pi*rand ;          % random initial position 
xy_=2*randn(2,1) ;  
xy=[[cos(phi_) -sin(phi_)] ; [sin(phi_) cos(phi_)]]*xy1+xy_ ; 
xy=xy+[Lx/2 ; Ly/2] ;     % adjust to frame center 
                          % ground truth emitter locaitons 

%%% Emitter activations parameters: Continuous illumination, state 0, 1, 2 
t0=0.40 ;                 % mean of T0 in sec
p0=exp(-Dt/t0)            % 0.9753, probability to retain state 0 in Dt
t1=0.025 ;                % mean of T1 in sec
p1=exp(-Dt/t1)            % 0.6703, probability to retain state 1 in Dt
t=4.5 ;                   % mean of T 
p=exp(-Dt/t)              % 0.9978, probability being photoactivatable in Dt
% transition probabilities
r00=p0*p                  % 0.9732, Trans. prob. from state 0 to 0
r10=(1-p0)*p              % 0.0246, Trans. prob. from state 0 to 1
r20=1-p                   % 0.0022, Trans. prob. from state 0 to 2
r01=(1-p1)*p              % 0.3290, Trans. prob. from state 1 to 0
r11=p1*p                  % 0.6688, Trans. prob. from state 1 to 1
r21=1-p                   % 0.0022, Trans. prob. from state 1 to 2

%%% Generate states of all emitters in data movie 
% and show number of activated emitters vs frame n
N=500 ;                   % number of frames in data movie
                          % Temporal resolution (TR): N*Dt=5 sec
[ca,p,h1]=emActMarkovContinue(t,t1,t0,Dt,N,M) ;
                          % h1=0.0697
a=(ca==1) ;               % a(m,n)=1 if activated; a(m,n)=0 otherwise  
                          % sum(sum(a')==0): # of emitters never activated 
Ma_=sum(a) ;              % number of activated emitters in nth frame
n=1:N ;
Ma=M*p.^n*h1 ;            % average # A-frames for M emitters in nth frame
% compare formulas and estimates 
Np=p*(1-p^N)/(1-p)        % 301.53
Np_=sum(sum(ca~=2))/M     % 302.03
Nae1=p*h1                 % 0.0695
Nae=((1-p^N)/(1-p))*Nae1  % 21.01
Nae_=sum(Ma_)/M           % 21.40
Na=M*Nae                  % 10504
Na_=M*Nae_                % 10699

%%% Fig. 11: number of activated emitters Ma(n) vs frame n
fprintf(1,'Fig. 11 \n') ;
wx=2.6 ; wy=0.75*wx ; dx=0.4 ; dy=0.34 ; 
PS=[dx/2+0.1 dy/2+0.1 wx-dx wy-dy] ;    % Positions of subfigures 
figure('Units','inches','Position',[3 4 wx wy],'Color',[1 1 1]) ;
subplot(1,1,1,'Units','inches','Position',PS) ;
plot(n,Ma_,'k.',n,Ma,'k-','MarkerSize',5) ;
axis([1 N 0 50])
ylabel('Number of activated emitters','FontSize',6)
xlabel('Frame n','FontSize',6) 
set(gca,'XTick',[0 50 100 150 200 250 300 350 400 450 500]) ;
set(gca,'XTickLabel',[0 50 100 150 200 250 300 350 400 450 500]) ;
set(gca,'YTick',[0 5 10 15 20 25 30 35 40 45 50]) ;
set(gca,'YTickLabel',[0 5 10 15 20 25 30 35 40 45 50]) ;
set(gca,'FontSize',6) ;
grid on
print('Fig11.tif','-dtiffn')

%%% Stack activated emitters in each frame
xyActive=zeros(2,M,N) ;  % activated emitter locations in nth frame
for n=1:N
  if Ma_(n)>0
    k=0 ; 
    for m=1:M
      if a(m,n)
        k=k+1 ;
        xyActive(:,k,n)=xy(:,m) ;
      end
      if k==Ma_(n), break ; end
    end
  end
end

%%% Fig. 12: Data frames 101, 102, 103, 104
fprintf(1,'Fig. 12 \n') ;
n12=101:104 ; 
em=[18 132 209 370 458] ;           % indices of emitters to be circled
st=['go' ; 'bo' ; 'mo' ; 'yo' ; 'co'] ; 
wx=1.3 ; wy=1.3 ; dx=0.04 ; dy=0.04 ; 
figure('Units','inches','Position',[1 2 4*wx wy],'Color',[1 1 1]) ;
PS=[0*wx+dx/2 dy/2 wx-dx wy-dy      % Positions of subfigures 
    1*wx+dx/2 dy/2 wx-dx wy-dy
    2*wx+dx/2 dy/2 wx-dx wy-dy
    3*wx+dx/2 dy/2 wx-dx wy-dy] ;
for i=1:length(n12)
  n=n12(i) ; 
  fprintf(1,'n=%d  Ma_(n)=%d \n',n,Ma_(n)) ; 
  k=Ma_(n) ;
  if k>0
    U=Gauss2D_Frame(sigma,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,xyActive(:,1:k,n)) ;
  else
    U=zeros(Ky,Kx) ;
  end
  subplot(1,4,i,'Units','inches','Position',PS(i,:)) ;
  show8bimage(U,'Yes','gray','No') ; hold on
  plot(xy(1,:)/Dx+0.5,xy(2,:)/Dy+0.5,'w.',xyActive(1,1:Ma_(n),n)/Dx+0.5, ... 
    xyActive(2,1:Ma_(n),n)/Dy+0.5,'r.','MarkerSize',3.5) ;
  for m=1:length(em)
    plot(xy(1,em(m))/Dx+0.5,xy(2,em(m))/Dy+0.5,st(m,:),'MarkerSize',3.5) ;
  end
  if i==1
      plot([2 12],Ky-[1 1],'w-',[2 2],Ky-[1-0.8 1+0.8],'w-',[12 12],Ky-[1-0.8 1+0.8],'w-') ;
      text(3.5,Ky-3.5,'1 {\mu}m','Color','white','FontSize',6) ;
  end
  hold off
  fm=compose('Frame %d',n) ; 
  text(27,3,fm,'Color','white','fontsize',6) ;
  axis off
  [min(min(U)) max(max(U))]
end
print('Fig12.tif','-dtiffn')
[[n12 ; a(:,n12)] [0 ; sum(a(:,n12),2)] (0:N)']   
Ma_(n12)   % [32    34    39    39]

%%% Generate a video of data movie: 100 frames  
fprintf(1,'Generate a video of data movie: 100 frames \n') ; 
vidObj=VideoWriter('Continuous2D.avi') ; 
vidObj.FrameRate=5 ;      % frames/s in movie
open(vidObj) ; 
wx=3.4 ; wy=wx ; dx=0.1 ; dy=0.1 ; 
PS=[dx/2 dy/2 wx-dx wy-dy] ;    % Positions of subfigures 
figure('Units','inches','Position',[3 4 wx wy],'Color',[1 1 1]) ;
subplot(1,1,1,'Units','inches','Position',PS) ;
for n=1:100
  k=Ma_(n) ;
  if k>0
    U=Gauss2D_Frame(sigma,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,xyActive(:,1:k,n)) ;
  else
    U=zeros(Ky,Kx) ;
  end
  show8bimage(U,'Yes','gray','No') ; hold on
  plot(xy(1,:)/Dx+0.5,xy(2,:)/Dy+0.5,'w.',xyActive(1,1:k,n)/Dx+0.5,xyActive(2,1:k,n)/Dy+0.5,'r.') ;
  plot([2 7],Ky-[1 1],'w-',[2 2],Ky-[1-0.5 1+0.5],'w-',[7 7],Ky-[1-0.5 1+0.5],'w-') ;
  hold off
  text(2.2,Ky-2.3,'500 nm','Color','white','FontSize',8) ;
  currFrame=getframe(gcf) ;
  writeVideo(vidObj,currFrame) ;
  fprintf(1,'N=%d n=%3d Ma_(n)=%3d \n',N,n,k) ;
end
close(vidObj) ;

%%% Generate and save a data movie
fprintf(1,'Generate a data movie: \n') ; 
for n=1:N
  k=Ma_(n) ;
  if k>0
    U=Gauss2D_Frame(sigma,Kx,Ky,Dx,Dy,Dt,Ih,b,mu,G,xyActive(:,1:k,n)) ;
  else
    U=zeros(Ky,Kx) ; 
  end
  U16=uint16(U) ;
  filename_Frame=strcat('Movie-Continuous2D\Continuous2D-',num2str(n),'.tif') ;
  imwrite(U16,filename_Frame) ;   % save data frame
  if mod(n,10)==0||n==1
    fprintf(1,'N=%3d n=%3d Ma_(n)=%d \n',N,n,Ma_(n)) ;
  end
end

%%% Verify saved frames 101-104
fprintf(1,'Verify saved frames 101-104 \n') ;
n12=101:104 ; 
em=[18 132 209 370 458] ;           % indices of emitters to be circled
st=['go' ; 'bo' ; 'mo' ; 'yo' ; 'co'] ; 
wx=1.7 ; wy=1.7 ; dx=0.04 ; dy=0.04 ; 
figure('Units','inches','Position',[1 2 4*wx wy],'Color',[1 1 1]) ;
PS=[0*wx+dx/2 dy/2 wx-dx wy-dy      % Positions of subfigures 
    1*wx+dx/2 dy/2 wx-dx wy-dy
    2*wx+dx/2 dy/2 wx-dx wy-dy
    3*wx+dx/2 dy/2 wx-dx wy-dy] ;
for i=1:length(n12)
  n=n12(i) ; 
  fprintf(1,'n=%d  Na(n)=%d \n',n,Ma_(n)) ; 
  filename_Frame=strcat('Movie-Continuous2D\Continuous2D-',num2str(n),'.tif') ;
  U16=imread(filename_Frame); % read a frame
  U=double(U16) ;
  subplot(1,4,i,'Units','inches','Position',PS(i,:)) ;
  show8bimage(U,'Yes','gray','No') ; hold on
  plot(xy(1,:)/Dx+0.5,xy(2,:)/Dy+0.5,'w.',xyActive(1,1:Ma_(n),n)/Dx+0.5, ...
    xyActive(2,1:Ma_(n),n)/Dy+0.5,'r.','MarkerSize',3.5) ;
  for m=1:length(em)
    plot(xy(1,em(m))/Dx+0.5,xy(2,em(m))/Dy+0.5,st(m,:),'MarkerSize',3.5) ;
  end
  if i==1
    plot([2 7],Ky-[1 1],'w-',[2 2],Ky-[1-0.5 1+0.5],'w-',[7 7],Ky-[1-0.5 1+0.5],'w-') ;
    text(1.8,Ky-2.5,'500 nm','Color','white','FontSize',5) ;
  end
  hold off
  fm=compose('Frame %d',n) ; 
  text(29.5,3,fm,'Color','white','fontsize',6) ;
  axis off
  [min(min(U)) max(max(U))]
end
[[n12 ; a(:,n12)] [0 ; sum(a(:,n12),2)] (0:N)']   
Ma_(n12)   % [32    34    39    39]
