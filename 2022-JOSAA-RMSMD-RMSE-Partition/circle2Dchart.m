% [xy,d]=circle2Dchart(R,M) 
%
% Generate a 2D circle chart for spatial resolution test.
% 
% Input:
%  R  - nm, resolution = distance between emitters
%  M  - number of emitters on circle
% Output:
%  xy	- size(2,N), xy(:,n) is nth vertex coordinate in x and y
%       dimension where N is total number of emitters in circle chart. 
%       The circle center is at [0 ; 0].
%  d	- diameter 
%
% 05/10/2020
% Yi Sun

function [xy,d]=circle2Dchart(R,M) 

%Dd=2*resol ; % nm, distance between adjacent emitters 
rd=R/(2*sin(pi/M)) ; % radius
theta=2*pi*(0:M-1)/M ; 
xy=[rd.*cos(theta) ; rd.*sin(theta)] ; 
d=2*rd ; 

end
