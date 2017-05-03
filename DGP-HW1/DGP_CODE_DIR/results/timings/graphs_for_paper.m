%timing graph:
MOSEK_timings=[764,3165,17517,101737,294947,723958,1645529,NaN,NaN,NaN];
MAP_gpu_timings=[77,87,188,256,379,666,882,1257,1567,2074];
ATP_timings=[17,23,40,72,132,225,318,436,562,760];
m=[500,1000,2000,4000,6000,8000,10000,12000,14000,16000];

figure
set(gca,'fontsize',12)
semilogy(m,ATP_timings,'LineWidth',4,'DisplayName','ATP');
hold on
semilogy(m,MAP_gpu_timings,'Color','r','LineWidth',4,'DisplayName','MAP-gpu');
hold on
semilogy(m,MOSEK_timings,'Color','c','LineWidth',4,'DisplayName','MOSEK');

legend('show')
% xlim([500 16000])
% ylim([0 2500])
xlabel('m')
ylabel('time (millisecond)')

%%
% 
%iterations graph:

MAP_gpu_timings=[203,176,199,173,141,143,135,133,131,130];
ATP_timings=[23,28,30,30,30,31,31,30,30,30];
m=[500,1000,2000,4000,6000,8000,10000,12000,14000,16000];

figure
set(gca,'fontsize',12)
semilogy(m,ATP_timings,'LineWidth',4,'DisplayName','ATP');
hold on
semilogy(m,MAP_gpu_timings,'Color','r','LineWidth',4,'DisplayName','MAP-gpu');

legend('show')
% xlim([500 16000])
% ylim([0 2500])
xlabel('m')
ylabel('Iterations')