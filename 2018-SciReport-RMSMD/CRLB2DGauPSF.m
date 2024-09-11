% [SDV SDVavg FWHM FWHMavg F_]=CRLB2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xy) 
%
% Calculate Cramer-Rao lower bound (CRLB) for 2-D nanoscopy imaging with a Gaussian 
% PSF. 
% 
% For technical details, see:
% [1] Sun, Y. Localization precision of stochastic optical localization
% nanoscopy using single frames. J. Biomed. Optics 18, 111418.1-15 (2013) 
% 
% Input:
%   sigma   - SDV of Gaussian PSF 
%		Kx, Ky	- Image size is Ky*Kx in pixels, sample is located at [0,Kx*Dx]x[0,Ky*Dy]
%		Dx, Dy	- Pixel size is Dx*Dy square nanometers (nm^2)
%   Dt   		- Imaging time in second (frame rate = 1/Dt)
%   Ih		  - Average number of detected photons per second per emitter 
%   SPNR 		- SPNR=1e-6*rp (um^2/emitter); 
%   SGNR 		- SGNR=1e-6*rg (um^2/emitter); 
%							Both SNRs can be equal to "inf"
%   xy      - ith colume is 2D coordinate (x,y)' of ith emitter location;
%
% Output:
%   SDV     - ith column of SDV is the standard deviation of estimate of ith column of xy
% 	SDVavg	- Average SDV of (x,y)
% 	FWHM		- FWHM=sqrt(2*log(2))*SDV 
% 	FWHMavg	- FWHMavg=sqrt(2*log(2))*SDVavg 
%   F_      - inverse of Fisher information matrix 
%
% Note:       All distances are in nm 
%             Gaussian noise is approximated by Poisson noise
%
% Yi Sun
% 1/23/2013, 3/27/2013, 7/27/2013

function [SDV SDVavg FWHM FWHMavg F_]=CRLB2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xy) 

[tmp M]=size(xy) ;
if M<1, 
  fprintf(1,'# of emitters is zero. \n') ;
  return ;
end
SNR=1/(1/SPNR+1/SGNR) ; r=10^6*SNR ; 
Q=zeros(Ky,Kx) ; 
qx=zeros(M,Kx) ;    qy=zeros(M,Ky) ; 
DqyDy=zeros(M,Ky) ; DqxDx=zeros(M,Kx) ; 
kx=0:Kx-1 ; ky=0:Ky-1 ; 
for i=1:M,
	Dxk1=(Dx*(kx+1)-xy(1,i))/sigma ; Dxk0=(Dx*kx-xy(1,i))/sigma ; 
	Dyk1=(Dy*(ky+1)-xy(2,i))/sigma ; Dyk0=(Dy*ky-xy(2,i))/sigma ; 
	qx(i,:)=(Qfunc(-Dxk1)-Qfunc(-Dxk0))/Dx ;
	qy(i,:)=(Qfunc(-Dyk1)-Qfunc(-Dyk0))/Dy ;
	Q=Q+qy(i,:)'*qx(i,:) ;
  DqyDy(i,:)=(1/Dy)*(-exp(-Dyk1.^2/2)+exp(-Dyk0.^2/2))/(sqrt(2*pi)*sigma) ;
  DqxDx(i,:)=(1/Dx)*(-exp(-Dxk1.^2/2)+exp(-Dxk0.^2/2))/(sqrt(2*pi)*sigma) ;
end
Qu=Q+1/r ;
if r==inf, 
  Qu=Qu+1e-10 ; % add 1e-10 to avoid 0/0=NaN 
end
F=zeros(2*M,2*M) ;	% Fisher information matrix of x1, y1, x2, y2, ..., xM, yM 
for i=1:M,	% 
	qyDqxDxi=qy(i,:)'*DqxDx(i,:) ; 
	DqyDyqxi=DqyDy(i,:)'*qx(i,:) ;
	for j=i:M,
		qyDqxDxj=qy(j,:)'*DqxDx(j,:) ;
		DqyDyqxj=DqyDy(j,:)'*qx(j,:) ;
% F(xi,:)
		F(2*(i-1)+1,2*(j-1)+1)=sum(sum(qyDqxDxi.*qyDqxDxj./Qu)) ;	% F(xi, xj)
		F(2*(j-1)+1,2*(i-1)+1)=F(2*(i-1)+1,2*(j-1)+1) ;						% F(xj, xi)
		F(2*(i-1)+1,2*(j-1)+2)=sum(sum(qyDqxDxi.*DqyDyqxj./Qu)) ;	% F(xi, yj)
		F(2*(j-1)+2,2*(i-1)+1)=F(2*(i-1)+1,2*(j-1)+2) ;						% F(yj, xi)
% F(yi,:)
		F(2*(i-1)+2,2*(j-1)+1)=sum(sum(DqyDyqxi.*qyDqxDxj./Qu)) ;	% F(yi, xj)
		F(2*(j-1)+1,2*(i-1)+2)=F(2*(i-1)+2,2*(j-1)+1) ;						% F(xj, yi)
		F(2*(i-1)+2,2*(j-1)+2)=sum(sum(DqyDyqxi.*DqyDyqxj./Qu)) ;	% F(yi, yj)
		F(2*(j-1)+2,2*(i-1)+2)=F(2*(i-1)+2,2*(j-1)+2) ;						% F(yj, yi)
	end
end
F=Ih*Dt*Dx*Dy*F ; F_=inv(F) ; 
var=diag(F_) ;
SDV=zeros(2,M) ;
for i=1:M,
	SDV(1,i)=sqrt(var(2*(i-1)+1)) ;	% SDV for x
	SDV(2,i)=sqrt(var(2*(i-1)+2)) ;	% SDV for y	
end
SDVavg=zeros(2,1) ;
SDVavg(1)=sum(SDV(1,:))/M ;
SDVavg(2)=sum(SDV(2,:))/M ;
FWHM=2*sqrt(2*log(2))*SDV ;
FWHMavg=2*sqrt(2*log(2))*SDVavg ;

end
