load("min_master_list_00.mat");
min_master_list_00 = min_master_list;

load("min_master_list_0.mat");
min_master_list_0 = min_master_list;

load("min_master_list_4.mat");
min_master_list_4 = min_master_list;

load("fine_min_master_list_8.mat");
min_master_list_8 = min_master_list;

load("fine_min_master_list_12.mat");
min_master_list_12 = min_master_list; 
load("fine_two_min_master_list_12.mat");
min_master_list_12 = cat(1,min_master_list_12,min_master_list);

max_power_list = -25:5:35;
plot(max_power_list,mean(min_master_list_12(:,1:end-1)),'DisplayName','12 IRS elements','Marker','x');
hold on;
% plot(max_power_list,mean(min_master_list_00),'DisplayName','No IRS with maximum power','Marker','>');
plot(max_power_list,mean(min_master_list_4(:,1:end-1)),'DisplayName','4 IRS elements','Marker','*');
% plot(max_power_list,mean(min_master_list_8),'DisplayName','8 IRS elements','Marker','^');
plot(max_power_list,mean(min_master_list_0(:,1:end-1)),'DisplayName','No IRS with optimized power','Marker','o');

ax = gca;
ax.FontSize = 12;
xlim([-25 35])
xlabel('Maximum transmit power (dBm)')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')
lgd = legend('Location','eastoutside');
lgd.FontSize = 12;
