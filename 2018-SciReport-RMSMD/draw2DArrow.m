% draw2DArrow(X0,X1,arrowSize,lineType,LineWidth) 
%
% Draw 2-D arrows in direction from X0 to X1 at X1.
% 
% Input:
%  X0        - size(X0)=[2,N]; 
%  X1        - size(X1)=[2,N]; an arrow is drawn for each pair from X0(:,i) to 
%              X1(:,i) at X1(:,i)
%  arrowSize - arrow size, arrowSize=[Size1, Size2]
%  arrowAngle - arrow angle in pi
%  lineType  - string. Line type is same as plot
%  LineWidth - default is 0.5 (point)
%
% 04/12/2017
% Yi Sun

function draw2DArrow(X0,X1,arrowSize,arrowAngle,lineType,LineWidth) 

[dim N0]=size(X0) ; % dim - dimension of X0, N0 - # of points in X0
[din N1]=size(X1) ; % din - dimension of X1, N1 - # of points in X1
if dim~=din, 
  fprintf(1,'X0 and X1 have different dimensions!\n')
  return ;
end
if N0~=N1,
  fprintf(1,'X0 and X1 have different columns!\n') ;
  return ;
end
if N0==0,
  fprintf(1,'X0 and X1 are empty!\n') ;
  return ;
end
angle1=pi-arrowAngle ; angle2=pi+arrowAngle ; 
A=zeros(2,4) ; 
hold on
for i=1:N0,
  x0=X0(:,i) ; x1=X1(:,i) ; 
  theta=atan((x1(2)-x0(2))/(x1(1)-x0(1))) ;
  if x1(2)-x0(2)<0&&x1(1)-x0(1)<0,
    theta=theta+pi ;
  end
  if x1(2)-x0(2)>0&&x1(1)-x0(1)<0,
    theta=theta+pi ;
  end
  A(:,2)=arrowSize(1)*[cos(theta+angle1) ; sin(theta+angle1)] ; 
  A(:,4)=arrowSize(2)*[cos(theta+angle2) ; sin(theta+angle2)] ; 
  plot(A(1,1:2)+x1(1),A(2,1:2)+x1(2),lineType) ; 
  plot(A(1,3:4)+x1(1),A(2,3:4)+x1(2),lineType) ; 
  plot([X0(1,i) X1(1,i)],[X0(2,i) X1(2,i)],lineType,'LineWidth',LineWidth) ;
end

end
