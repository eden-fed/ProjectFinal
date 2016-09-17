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
x = 0:0.01:0.45;
y = -log((1-abs(x))/0.5);
figure;
plot(x,y,'LineWidth',3)

hold on

x = 0:0.01:0.45;
y = 0.7-x;
plot(x,y,'LineWidth',3)
line([0.4 0.4],[-1 1.5],'LineWidth',3,'Color','m');



%**********************************************************



plot(abs_Vg1,R_l_gz1,'LineWidth',3,'g');
    hold on
    plot(abs_Vg,R_l_gz,'LineWidth',2,'y');
    hold on ;
		x = 0:0.01:0.45;
y = -log((1-abs(x))/0.5);
plot(x,y,'LineWidth',3)

hold on

x = 0:0.01:0.45;
y = 0.7-x;
plot(x,y,'LineWidth',3)
line([0.4 0.4],[-1 1.5],'LineWidth',3,'Color','m');