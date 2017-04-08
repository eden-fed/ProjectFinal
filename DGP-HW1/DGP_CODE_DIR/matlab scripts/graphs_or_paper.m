%timing graph:
cvx_timings=[764,3165,17517,101737,294947,NaN,NaN,NaN,NaN,NaN];
LG_gpu_timings=[77,87,188,256,379,666,882,1257,1567,2074];
ATP_timings=[17,23,40,72,132,225,318,436,562,760];
m=[500,1000,2000,4000,6000,8000,10000,12000,14000,16000];

figure
set(gca,'fontsize',12)
plot(m,ATP_timings,'LineWidth',4,'DisplayName','ATP');
hold on
plot(m,LG_gpu_timings,'Color','r','LineWidth',4,'DisplayName','LG gpu');
hold on
plot(m,cvx_timings,'Color','c','LineWidth',4,'DisplayName','CVX');

legend('show')
xlim([500 16000])
ylim([0 2500])
xlabel('m (n=0.1*m)')
ylabel('time (millisecond)')

%%

%iterations graph:

LG_gpu_timings=[203,176,199,173,141,143,135,133,131,130];
ATP_timings=[23,28,30,30,30,31,31,30,30,30];
m=[500,1000,2000,4000,6000,8000,10000,12000,14000,16000];

figure
set(gca,'fontsize',12)
plot(m,ATP_timings,'LineWidth',4,'DisplayName','ATP');
hold on
plot(m,LG_gpu_timings,'Color','r','LineWidth',4,'DisplayName','LG gpu');

legend('show')
% xlim([500 16000])
% ylim([0 2500])
xlabel('m (n=0.1*m)')
ylabel('time (millisecond)')