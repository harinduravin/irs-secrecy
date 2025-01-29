load("min_master_list.mat");

load("all_master_list.mat");

max_power_list = -25:5:40;

mean_user_rates = mean(all_master_list,3);

figure(1)
plot(max_power_list,mean(min_master_list,'omitnan'),'DisplayName','Minimum secrecy rate','Marker','x');
hold on;
plot(max_power_list,mean_user_rates(:,1)','DisplayName','A1','Marker','o');
plot(max_power_list,mean_user_rates(:,2)','DisplayName','B1','Marker','o');
plot(max_power_list,mean_user_rates(:,3)','DisplayName','A2','Marker','o');
plot(max_power_list,mean_user_rates(:,4)','DisplayName','B2','Marker','o');
plot(max_power_list,mean_user_rates(:,5)','DisplayName','A3','Marker','o');
plot(max_power_list,mean_user_rates(:,6)','DisplayName','B3','Marker','o');

ax = gca;
ax.FontSize = 12;
xlim([-25 35])
xlabel('Maximum transmit power (dBm)')
ylabel('Secrecy Rate(bits/sec/Hz)')
title('Maximum transmit power against secrecy rate (4 IRS elements)')
lgd = legend('Location','eastoutside');
lgd.FontSize = 12;


