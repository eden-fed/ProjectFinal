x = 0:0.01:1;
y = -log((1-abs(x))/0.5);
figure;
plot(x,y,'LineWidth',3);

hold on;


x = 0:0.124:1;
y = -log((1-abs(x))/0.5);
plot(x,y,'-o','LineWidth',2);

%**********************************************************
x = 0:0.01:0.8;
y = -log((1-abs(x))/0.5);
figure;
plot(x,y,'LineWidth',3)

hold on

x = 0:0.01:0.8;
y = 0.7-x;
plot(x,y,'LineWidth',3)
line([0.4 0.4],[-1 1.5],'LineWidth',3,'Color','m');

%**********************************************************
%draw the region:
sigma=0.5;
SIGMA=2;
k=0.4;

x = 0:0.01:0.7;
y = -log((1-x)/sigma);
figure;
plot(x,y,'LineWidth',3)

hold on

y = log(SIGMA)-x;
plot(x,y,'LineWidth',3)
line([k k],[-1 1.5],'LineWidth',3,'Color','m');
line([0 0.7],[0 0],'LineWidth',1,'Color','k');

axis equal tight 

%draw the region with non-convex part:
% sigma=0.5;
% SIGMA=2;
% k=0.4;

x = 0:0.01:0.7;
y = -log((1-x)/sigma);
figure;
plot(x,y,'LineWidth',3)

hold on
y=-log((1+x)/SIGMA);
plot(x,y,'LineWidth',3,'color','g');
hold on

y = log(SIGMA)-x;
plot(x,y,'LineWidth',3)
line([k k],[-1 1.5],'LineWidth',3,'Color','c');
line([0 0.7],[0 0],'LineWidth',1,'Color','k');

% hold on
% plot(abs(Vg),real(l_gz),'.','Color','r')
% xlabel('abs(V_gz)');
% ylabel('real(L_gz)');

axis equal tight 
%**********************************************************



plot(abs_Vg1,R_l_gz1,'LineWidth',3,'Color','g');
hold on
plot(abs_Vg,R_l_gz,'LineWidth',2,'Color','b');
hold on ;
	x = 0:0.01:0.45;
y = -log((1-abs(x))/0.5);
plot(x,y,'LineWidth',3,'Color','m')

hold on

x = 0:0.01:0.45;
y = 0.7-x;
plot(x,y,'LineWidth',3,'Color','m')
line([0.4 0.4],[-1 1.5],'LineWidth',3,'Color','m');


%*************************************************************


x = 0:0.01:0.45;
y = -log((1-abs(x))/0.5);
plot(x,y,'LineWidth',3,'Color','m')

hold on

x = 0:0.01:0.45;
y = 0.7-x;
plot(x,y,'LineWidth',3,'Color','m')
line([0.4 0.4],[-1 1.5],'LineWidth',3,'Color','m');
hold on
plot(abs_Vg1(234),R_l_gz1(234),'g*');
hold on
plot(abs_Vg(234),R_l_gz(234),'b*');
hold on ;
plot(inputSampleX,inputSampleY,'c*')

%********************************************************
temp=(abs_Vg1-abs_Vg).^2+(R_l_gz1-R_l_gz).^2
[~,I]=max(temp)

%*********************************************************
%draw the region with non-convex part:


x = -0.9:0.01:0.9;
y = -log((1-abs(x))/sigma);
figure;
plot(x,y,'LineWidth',3,'color','g')

hold on
y=-log((1+abs(x))/SIGMA);
plot(x,y,'LineWidth',3,'color','g');
hold on

y = log(SIGMA)-abs(x);
plot(x,y,'LineWidth',3)
line([k k],[-0.7 1.7],'LineWidth',3,'Color','c');
line([-k -k],[-0.7 1.7],'LineWidth',3,'Color','c');
line([-0.9 0.9],[0 0],'LineWidth',1,'Color','k');

hold on

y=m*abs(x)+b;
plot(x,y,'LineWidth',3);

axis equal tight 
%******************************************************
%draw the region:
sigma=0.5;
SIGMA=2;
k=0.4;

x = 0:0.01:0.7;
figure;

hold on

y = log(SIGMA)-x;
plot(x,y,'LineWidth',3)
line([k k],[-1 1.5],'LineWidth',3,'Color','m');
line([0 0.7],[0 0],'LineWidth',1,'Color','k');

hold on

y=m*abs(x)+b;
plot(x,y,'LineWidth',3);

axis equal tight 

%***********************************************************
sigma=0.5;
SIGMA=2;
k=0.4;
crossPointY=log(sigma/(1-k));
m=(crossPointY-log(sigma))/k;
b=log(sigma);

x = 0:0.01:k;
figure;

hold on

y = log(SIGMA)-x;
plot(x,y,'LineWidth',3)
line([k k],[log(sigma)-0.3 log(SIGMA)+0.4],'LineWidth',3);
line([0 0.7],[0 0],'LineWidth',1,'Color','k');

hold on

y=m*x+b;
plot(x,y,'LineWidth',3);

hold on

y=x+log(SIGMA);
plot(x,y,'LineWidth',3,'Color','m');

hold on

y=-x/m+log(sigma);
plot(x,y,'LineWidth',3,'Color','m');

hold on

x = k:0.01:2*k;
y=-x/m+crossPointY+k/m;
plot(x,y,'LineWidth',3,'Color','m');

hold on

y=x+log(SIGMA)-2*k;
plot(x,y,'LineWidth',3,'Color','m');

hold on

y=log(SIGMA)-k;
plot(x,y,'LineWidth',3,'Color','m');

hold on

y=log(sigma/(1-k));
plot(x,y,'LineWidth',3,'Color','m');

axis equal tight 
%***************************************************************
sigma=0.5;
SIGMA=2;
k=0.7;
crossPointX=-(lambertw((sigma/SIGMA)*exp(1))-1);
crossPointY=log(sigma/(1-crossPointX));
m=(crossPointY-log(sigma))/crossPointX;
b=log(sigma);

x = 0:0.01:2*crossPointX;
figure;

hold on

y = log(SIGMA)-x;
plot(x,y,'LineWidth',3)
line([k k],[log(sigma)-0.3 log(SIGMA)+0.4],'LineWidth',1);
line([0 0.7],[0 0],'LineWidth',1,'Color','k');

hold on

y=m*x+b;
plot(x,y,'LineWidth',3);

hold on

x = 0:0.01:crossPointX;
y=x+log(SIGMA);
plot(x,y,'LineWidth',3,'Color','m');

hold on

y=-x/m+log(sigma);
plot(x,y,'LineWidth',3,'Color','m');

hold on

x = crossPointX:0.01:2*crossPointX;
y=x+crossPointY-crossPointX;
plot(x,y,'LineWidth',3,'Color','m');

hold on

y=-x/m+crossPointY+crossPointX/m;
plot(x,y,'LineWidth',3,'Color','m');

axis equal tight 