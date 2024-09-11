% Q-function, complementary cumulatve 
% Gaussian distribution function
% y=0.5*(1+erf(-x/2^0.5)) ;
%
% Yi Sun
% 12/8/2005

function y=Qfunc(x)
 y=0.5*(1+erf(-x/2^0.5)) ;
 z=(x>7) ;
 if sum(z)==0,
  return ;
 else
  for i=1:length(x),
   if z(i)==1,
    y(i)=exp(-x(i)^2/2)/(2*pi)^0.5/x(i)*(1-1/x(i)^2+3/x(i)^4) ;
   end
  end
 end






