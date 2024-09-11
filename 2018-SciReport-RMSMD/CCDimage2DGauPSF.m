% V=CCDimage2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xy) 
%
% Produce a CCD image (dada frame) with a 2D Gaussian PSF 
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
%             xy can be empty, i.e. xy=[]
%
% Output:
%   V       - a CCD image (data frame)
%
% Note:       Assume mean of Gaussian noise is zero
%             All distances are in nm
%
% Yi Sun
% 1/23/2013, 5/1/2018 

function V=CCDimage2DGauPSF(sigma,Kx,Ky,Dx,Dy,Dt,Ih,SPNR,SGNR,xy) 

[tmp M]=size(xy) ; 
rp=10^6*SPNR ; rg=10^6*SGNR ; G=Ih/rg ; 
Q=zeros(Ky,Kx) ; 
if M>=1,
  qx=zeros(M,Kx) ; qy=zeros(M,Ky) ;
  kx=0:Kx-1 ; ky=0:Ky-1 ; 
  for i=1:M,
    Dxk1=(Dx*(kx+1)-xy(1,i))/sigma ; Dxk0=(Dx*kx-xy(1,i))/sigma ;
    Dyk1=(Dy*(ky+1)-xy(2,i))/sigma ; Dyk0=(Dy*ky-xy(2,i))/sigma ;
    qx(i,:)=(Qfunc(-Dxk1)-Qfunc(-Dxk0)) ;
    qy(i,:)=(Qfunc(-Dyk1)-Qfunc(-Dyk0)) ;
    Q=Q+qy(i,:)'*qx(i,:) ;
  end
  Q=(Dx*Dy)^(-1)*Q ;
end
v=Ih*Dt*Dx*Dy*(Q+1/rp) ;
V=poissrnd(v) ;
V=V+sqrt(Dt*Dx*Dy*G)*randn(size(V)) ;

end
