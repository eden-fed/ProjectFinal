function [ h ] = convex_graph_with_map( sigma, SIGMA, k, l_gz, Vg ,energy,m,b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = 0:0.01:0.7;
% y = -log((1-x)/sigma);
% plot(x,y,'LineWidth',3,'color','g')

% hold on
% y=-log((1+x)/SIGMA);
% plot(x,y,'LineWidth',3,'color','g');
hold on

y = log(SIGMA)-x;
plot(x,y,'LineWidth',3,'color','b')

line([k k],[-1 1],'LineWidth',3,'Color','b');
line([0 0.7],[0 0],'LineWidth',1,'Color','k');

y=m*abs(x)+b;
plot(x,y,'LineWidth',3,'color','b');

hold on
h=plot(abs(Vg),real(l_gz),'.','Color','r','MarkerSize',20);
xlabel('abs(V_gz)');
ylabel('real(L_gz)');
title(['energy = ',num2str(energy)]);

axis equal tight 
drawnow

end

