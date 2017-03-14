function [ x_poly,y_poly ] = projectToPolywithK( x,y,SIGMA,sigma,k )
%x = vector of x coordinates of the given points
%y = vector of y coordinate of the given points

log_SIGMA=log(SIGMA);
log_sigma=log(sigma);

if(k==0)
    m=1;
else
    m=(log(sigma/(1-k))-log_sigma)/k;
end

%divide all the points to regions (according to picture of case 1)
I1 = (y>=x+log_SIGMA);
I2 = (x > log_SIGMA-y) & (y<x+log_SIGMA) & (y>x+log_SIGMA-2*k);
I3 = (y<=x+log_SIGMA-2*k) & (y>=log_SIGMA-k);
I4 = (x>k) & (y<log_SIGMA-k) & (y>log(sigma/(1-k)));
I5 = (y<=log(sigma/(1-k))) & (y>=-x/m+log(sigma/(1-k))+k/m);
I6 = (y<-x/m+log(sigma/(1-k))+k/m) & (y<m*x+log_sigma) & (y>-x/m+log_sigma);
I7 = (y<=-x/m+log_sigma);

x_poly=x;
y_poly=y;

%find the projection on the lines, for all the regions
x_poly(I1)=0;
y_poly(I1)=log_SIGMA;

x_poly(I2)=(x(I2)-y(I2)+log_SIGMA)/2;
y_poly(I2)=(y(I2)-x(I2)+log_SIGMA)/2;

x_poly(I3)=k;
y_poly(I3)=log_SIGMA-k;

x_poly(I4)=k;
y_poly(I4)=y(I4);

x_poly(I5)=k;
y_poly(I5)=log(sigma/(1-k));

x_poly(I6)=(x(I6)+m*y(I6)-m*log_sigma)/(m^2+1);
y_poly(I6)=m*((x(I6)+m*y(I6)-m*log_sigma)/(m^2+1))+log_sigma;

x_poly(I7)=0;
y_poly(I7)=log_sigma;

end