% [ca,p,h1]=emActMarkovCycle(t,tA1,tA0,tD1,tD0,Dt,C,KA,KD,M)
%
%   Produce photoactivation states for emitters in a data movie based on 
% two-phase Markov chain for cycled illumination. 
%   Each cycle consists of two phases: A-phase consists of KA frames and 
% D-phase consists of KD frames. 
% 
% Input:
%   t       - mean of photoactivatable period 
%   tA1     - mean of activation period in A-phase, tA1 << t
%   tA0     - mean of deactivation period in D-phase, tA0 << t, tA0 > tA1
%   tD1     - mean of activation period in D-phase, tD1 << t
%   tD0     - mean of deactivation period in D-phase, tD0 << t, tD0 < tD1
%           - tA0 << tD0, tA1 >> tD1
%   Dt      - frame time, Dt << tA0, tD0, t
%   C       - number of cycles
%   KA      - number of frames in A-phase of a cycle
%   KD      - number of frames in D-phase of a cycle
%           - K = KA + KD, number of frames per cycle
%           - N = C*K. number of frames in a data movie
%   M       - number of emitters 
% Output:
%   ca      - ca(M,N), state of mth emitter in nth frame ca(m,n): 
%             0 - deactivated
%             1 - activated
%             2 - photobleached
%   p       - probability that a frame is photoactivatable if it is at
%             beginning
%   h1      - h1(1,K), h1(i) is stationary probability of state 1 in ith
%             frame of a cycle as t -> inf
%
% Note:   All times are in second
%
% 04/26/2024

function [ca,p,h1]=emActMarkovCycle(t,tA1,tA0,tD1,tD0,Dt,C,KA,KD,M) 

K=KA+KD ;
N=C*K ; 
% A-phase
pA0=exp(-Dt/tA0) ;  % probability to retain state 0 in Dt
pA1=exp(-Dt/tA1) ;  % probability to retain state 1 in Dt
% D-phase 
pD0=exp(-Dt/tD0) ;  % probability to retain state 0 in Dt
pD1=exp(-Dt/tD1) ;  % probability to retain state 1 in Dt
% Bleaching
p=exp(-Dt/t) ;      % probability to retain photoactivatable in Dt
% transistion matrices as t -> inf
A=[pA0 1-pA1 ; 1-pA0 pA1] ; 
D=[pD0 1-pD1 ; 1-pD0 pD1] ; 
h1=zeros(1,K) ; 
for i=1:K
  if i<=KA
    Q=A^i*D^KD*A^(KA-i) ; 
  else
    Q=D^(i-KA)*A^KA*D^(K-i) ; 
  end
  h1(i)=(1-Q(1,1))/(2-Q(1,1)-Q(2,2)) ;
end
% transition probabilities
a00=pA0*p ;         
a10=(1-pA0)*p ;     
a20=1-p ;           
a01=(1-pA1)*p ;     
a11=pA1*p ;         
a21=1-p ;           
d00=pD0*p ;         
d10=(1-pD0)*p ;     
d20=1-p ;           
d01=(1-pD1)*p ;     
d11=pD1*p ;         
d21=1-p ;           
Ra=[a00 a01 0       % State transition matrix for A-phase
    a10 a11 0
    a20 a21 1] ;
Rd=[d00 d01 0       % State transition matrix for D-phase
    d10 d11 0
    d20 d21 1] ;
c0=zeros(M,N+1) ;   % states of Markov chains in data movie
c0(:,1)=(rand(M,1)<=h1(K)) ;  % initial states: 0, 1 in stationary distr
n=0 ; 
for nc=0:C-1        % C cycles in a data movie
  % A-phase
  for i=1:KA
    n=n+1 ;
    for m=1:M
      if c0(m,n)==2 % photobleached 
        c0(m,n+1)=2 ;
      else          % c0(m,n)==0 or 1
        P=rand ;
        if c0(m,n)==0
          if P<Ra(1,1), c0(m,n+1)=0 ; end
          if P>=Ra(1,1)&&P<Ra(1,1)+Ra(2,1), c0(m,n+1)=1 ; end
          if P>=Ra(1,1)+Ra(2,1), c0(m,n+1)=2 ; end
        else        % c0(m,n)==1
          if P<Ra(1,2), c0(m,n+1)=0 ; end
          if P>=Ra(1,2)&&P<Ra(1,2)+Ra(2,2), c0(m,n+1)=1 ; end
          if P>=Ra(1,2)+Ra(2,2), c0(m,n+1)=2 ; end
        end
      end
    end
  end
  % D-phase 
  for i=1:KD
    n=n+1 ;
    for m=1:M
      if c0(m,n)==2
        c0(m,n+1)=2 ;
      else          % c0(m,n)==0 or 1
        P=rand ;
        if c0(m,n)==0
          if P<Rd(1,1), c0(m,n+1)=0 ; end
          if P>=Rd(1,1)&&P<Rd(1,1)+Rd(2,1), c0(m,n+1)=1 ; end
          if P>=Rd(1,1)+Rd(2,1), c0(m,n+1)=2 ; end
        else        % c0(m,n)==1
          if P<Rd(1,2), c0(m,n+1)=0 ; end
          if P>=Rd(1,2)&&P<Rd(1,2)+Rd(2,2), c0(m,n+1)=1 ; end
          if P>=Rd(1,2)+Rd(2,2), c0(m,n+1)=2 ; end
        end
      end
    end
  end
end
ca=c0(:,2:N+1) ;    % remove initial state

end
