
load("min_master_list_0.mat");
min_master_list_0 = min_master_list;

load("min_master_list_4.mat");
min_master_list_4 = min_master_list;


load("fine_min_master_list_12.mat");
min_master_list_12 = min_master_list; 
load("fine_two_min_master_list_12.mat");
min_master_list_12 = cat(1,min_master_list_12,min_master_list);

max_power_list = -25:5:25;

plot(max_power_list,mean(min_master_list_0(:,1:end-3)),'DisplayName','No IRS','Marker','o','LineWidth',1.5);
hold on;
plot(max_power_list,mean(min_master_list_4(:,1:end-3)),'DisplayName','4 IRS elements','Marker','*','LineWidth',1.5);

plot(max_power_list,mean(min_master_list_12(:,1:end-3)),'DisplayName','12 IRS elements','Marker','x','LineWidth',1.5);

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
xlim([-25 25])
xlabel('Maximum transmit power (dB)','Interpreter','latex','FontName','Times','FontSize',12)
ylabel('Minimum secrecy Rate(bits/sec/Hz)','Interpreter','latex','FontName','Times','FontSize',12)
lgd = legend('Location','best','FontName','Times');
lgd.FontSize = 12;
grid on;

print -depsc2 chart_power.eps