close all
sigma=0.5;
SIGMA=2;
k=0.2+rand(1,1);
crossPointX=-(lambertw((sigma/SIGMA)*exp(1))-1);

if(k<crossPointX)
    xv=[0;k;k;0];
    yv=[log(SIGMA);log(SIGMA)-k;log(sigma/(1-k));log(sigma)];
else
    xv=[0;crossPointX;0];
    yv=[log(SIGMA);log(sigma/(1-crossPointX));log(sigma)];
end

a=20000;
x_in=rand(a,1);
y_in=-1+2*rand(a,1);
[m,intersectionX]=findLineForTest(sigma,SIGMA,k);

% tic
% x_out1=x_in;
% y_out1=y_in; 
% I=(x_in > log(SIGMA)-y_in) | (x_in>k) | (x_in<0) | (y_in<m*x_in+b);
% [x_poly,y_poly] = p_poly_dist_vec(x_in(I), y_in(I), xv, yv);
% x_out1(I)=x_poly;
% y_out1(I)=y_poly;
% using_p_poly_dist=toc

% tic
% x_out=zeros(size(x_in));
% y_out=zeros(size(x_in));
% for ii=1:length(x_in)
%     [x_poly,y_poly] = p_poly_dist(x_in(ii), y_in(ii), xv, yv);
%     x_out(ii)=x_poly;
%     y_out(ii)=y_poly;
% end
% using_loop_p_poly_dist=toc
tic
[ x_out,y_out ] = projectToPoly( x_in,y_in,SIGMA,sigma,k,intersectionX );
using_projectToPoly=toc

tic
 cvx_begin 
    variable x_out1(a);
    variable y_out1(a);
    minimize sum_square_abs(x_in-x_out1)+sum_square_abs(y_in-y_out1);
    subject to
        x_out1>=0;
        x_out1<=k;
        x_out1<=log(SIGMA)-y_out1;
        m*x_out1+log(sigma)<=y_out1;
cvx_end 
using_cvx=toc
    
% tic
% x_out=zeros(size(x_in));
% y_out=zeros(size(x_in));
% for ii=1:length(x_in)
%     [x_poly,y_poly] = projectPointToPoly(x_in(ii), y_in(ii), SIGMA,sigma,k);
%     x_out(ii)=x_poly;
%     y_out(ii)=y_poly;
% end
% using_loop_projectToPoly=toc
% 
tic
x_out2=zeros(size(x_in));
y_out2=zeros(size(x_in));
for ii=1:length(x_in)
    [x_poly,y_poly] = projectPointToPolyMinSeg(x_in(ii), y_in(ii), xv, yv);
    x_out2(ii)=x_poly;
    y_out2(ii)=y_poly;
end
using_loop_projectPointToPolyMinSeg=toc
% 

if(any(abs(x_out1-x_out) > 1e-5) || any(abs(y_out1-y_out) > 1e-5))
% if(any(sqrt((x_out1-x_out).^2+(y_out1-y_out).^2)>5e-4))
    [M,I]=max((x_out1-x_out).^2+(y_out1-y_out).^2);
    warning('the results are different when x=%d , y=%d',x_in(I),y_in(I));
end
    
% draw the region:

 crossPointY=log(sigma/(1-intersectionX));
m=(crossPointY-log(sigma))/intersectionX;

x = 0:0.01:intersectionX;
figure;

hold on

y = log(SIGMA)-x;
plot(x,y,'LineWidth',3)
line([k k],[log(sigma)-0.3 log(SIGMA)+0.4],'LineWidth',3);
line([0 0.7],[0 0],'LineWidth',1,'Color','k');

hold on

y=m*x+log(sigma);
plot(x,y,'LineWidth',3);

% hold on
% x = crossPointX:0.01:2.5*crossPointX;
% y=-x/m+crossPointY+crossPointX/m;
% plot(x,y,'LineWidth',3,'Color','m');
% 
% hold on
% y=log(sigma/(1-k));
% plot(x,y,'LineWidth',3,'Color','m');

hold on
plot(x_in,y_in,'.','Color','r')
hold on
plot(x_out,y_out,'.','Color','g')
hold on
plot(x_out1,y_out1,'.','Color','c')

axis equal tight 