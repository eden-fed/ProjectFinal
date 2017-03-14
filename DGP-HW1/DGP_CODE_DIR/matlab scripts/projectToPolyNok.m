function [ x_poly,y_poly ] = projectToPolyNok( x,y,SIGMA,sigma,intersectionX)
%x = vector of x coordinates of the given points
%y = vector of y coordinate of the given points
%intersectionX = x coordinate of the intersection point of the first and third constrains

log_SIGMA=log(SIGMA);
log_sigma=log(sigma);
intersectionY=log(sigma/(1-intersectionX));

if(intersectionX==0)
    m=1;
else
    m=(intersectionY-log_sigma)/intersectionX;
end

%divide all the points to regions (according to picture of case 2)
I1 = (y>=x+log_SIGMA);
I2 = (x > log_SIGMA-y) & (y<x+log_SIGMA) & (y>x+intersectionY-intersectionX);
I3 = (y<=x+intersectionY-intersectionX) & (y>=-x/m+intersectionY+intersectionX/m);
I4 = (y<-x/m+intersectionY+intersectionX/m) & (y<m*x+log_sigma) & (y>-x/m+log_sigma);
I5 = (y<=-x/m+log_sigma);

x_poly=x;
y_poly=y;

%find the projection on the lines, for all the regions
x_poly(I1)=0;
y_poly(I1)=log_SIGMA;

x_poly(I2)=(x(I2)-y(I2)+log_SIGMA)/2;
y_poly(I2)=(y(I2)-x(I2)+log_SIGMA)/2;

x_poly(I3)=intersectionX;
y_poly(I3)=intersectionY;

x_poly(I4)=(x(I4)+m*y(I4)-m*log_sigma)/(m^2+1);
y_poly(I4)=m*((x(I4)+m*y(I4)-m*log_sigma)/(m^2+1))+log_sigma;

x_poly(I5)=0;
y_poly(I5)=log_sigma;

end
