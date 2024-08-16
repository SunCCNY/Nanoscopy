% [ca,p,h1]=emActMarkovContinue(t,t1,t0,Dt,N,M)
%
% Produce emitters' photoactivation states in a data movie based on 
% two-phase Markov chain for cycled illumination. 
% Each cycle consists two phase: A-phase consists of KA frames and 
% D-phase consists of KD frames. 
% 
% Input:
%   t       - mean of photoactivatable period 
%   t1      - mean of activation period, t1 << t
%   t0      - mean of deactivation period, t0 << t, t0 > t1
%   Dt      - frame time, Dt << t0, t
%   N       - number of frames in a data movie
%   M       - number of emitters 
% Output:
%   ca      - ca(M,N), state of mth emitter in nth frame ca(m,n): 
%             0 - deactivated
%             1 - activated
%             2 - photobleached
%   p       - probability that a frame is photoactivatable if it is at
%             beginning
%   h1      - stationary probability of state 1 as t -> inf
%
% Note:   All times are in second
%
% 04/29/2024

function [ca,p,h1]=emActMarkovContinue(t,t1,t0,Dt,N,M)

p0=exp(-Dt/t0) ;    % probability to retain state 0 in Dt
p1=exp(-Dt/t1) ;    % probability to retain state 1 in Dt
p=exp(-Dt/t) ;      % probability to retain photoactivatable in Dt
h1=(1-p0)/(2-p0-p1) ;
% transition probabilities
r00=p0*p ;         
r10=(1-p0)*p ;     
r20=1-p ;           
r01=(1-p1)*p ;     
r11=p1*p ;         
r21=1-p ;           
R=[r00 r01 0        % State transition matrix 
   r10 r11 0
   r20 r21 1] ;
c0=zeros(M,N+1) ;   % states of Markov chains in data movie
c0(:,1)=(rand(M,1)<=h1) ;  % initial states: 0, 1 in stationary distr
for n=1:N
  for m=1:M
    if c0(m,n)==2   % photobleached 
      c0(m,n+1)=2 ; 
    else            % c0(m,n)==0 or 1
      P=rand ; 
      if c0(m,n)==0
        if P<R(1,1), c0(m,n+1)=0 ; end
        if P>=R(1,1)&&P<R(1,1)+R(2,1), c0(m,n+1)=1 ; end
        if P>=R(1,1)+R(2,1), c0(m,n+1)=2 ; end
      else          % c0(m,n)==1
        if P<R(1,2), c0(m,n+1)=0 ; end
        if P>=R(1,2)&&P<R(1,2)+R(2,2), c0(m,n+1)=1 ; end
        if P>=R(1,2)+R(2,2), c0(m,n+1)=2 ; end
      end        
    end
  end
end
ca=c0(:,2:N+1) ;   % remove initial all-0 state

end
