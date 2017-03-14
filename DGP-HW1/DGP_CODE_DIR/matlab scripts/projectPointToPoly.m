function [ x_poly,y_poly ] = projectPointToPoly( x,y,SIGMA,sigma,k )

if (x <= log(SIGMA)-y) && (x<=k) && (y>=log(sigma)-(log(1-k)/k)*x)
    x_poly=x;
    y_poly=y;
elseif (y>=x+log(SIGMA))
    x_poly=0;
    y_poly=log(SIGMA);
elseif (x > log(SIGMA)-y) && (y<x+log(SIGMA)) && (y>x+log(SIGMA)-2*k)
    x_poly=(x-y+log(SIGMA))/2;
    y_poly=(y-x+log(SIGMA))/2;
elseif (y<=x+log(SIGMA)-2*k) && (y>=log(SIGMA)-k)
    x_poly=k;
    y_poly=log(SIGMA)-k;
elseif (x>k) && (y<log(SIGMA)-k) && (y>log(sigma/(1-k)))
    x_poly=k;
    y_poly=y;
elseif (y<=log(sigma/(1-k))) && (y>=x*k/log(1-k)+log(sigma/(1-k))-(k^2)/log(1-k))
    x_poly=k;
    y_poly=log(sigma/(1-k));
elseif (y<=x*k/log(1-k)+log(sigma))
    x_poly=0;
    y_poly=-log(1/sigma);
else
    x_poly=(x*k-log(1-k)*y+log(1-k)*log(sigma))*k/((log(1-k))^2+k^2);
    y_poly=(log(sigma)*k^2+log(1-k)^2*y-log(1-k)*k*x)/((log(1-k))^2+k^2);
end

end