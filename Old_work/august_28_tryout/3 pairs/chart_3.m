load("3_pairs_20_seeds_2_trial.mat");
all_min_inf_master_list = min_inf_master_list;
all_min_leaked_master_list = min_leaked_master_list;
all_min_master_list = min_master_list;

load("3_pairs_20_seeds_1_trial.mat");
all_min_inf_master_list = cat(1,all_min_inf_master_list,min_inf_master_list);
all_min_leaked_master_list = cat(1,all_min_leaked_master_list,min_leaked_master_list);
all_min_master_list = cat(1,all_min_master_list,min_master_list);

mean_inf = arrayfun(@watt2dbm, mean(all_min_inf_master_list));

L_list = 0:5:40;
L_list(1) = 1;

figure (2)
yyaxis left
plot(L_list,max(mean(all_min_master_list),0),'DisplayName','Minimum Secrecy Rate','LineWidth',1.5,'Marker','^','Color',[0.9 0.6 0]);
hold on;
plot(L_list,mean(all_min_leaked_master_list),'DisplayName','Leaked Information Rate','LineWidth',1.5);
xlabel('Number of IRS elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')
lgd = legend('Location','best');

yyaxis right
plot(L_list,mean_inf,'DisplayName','Interference power','LineWidth',1.5);
ylabel('Experienced Interference Power (dBm)')
title({'Interference and leaked information rate of','the most vulnerable user (3 user pairs in the system)'})

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
lgd.FontSize = 12;
lgd.Position = [0.499180963622548,0.158422620413559,0.378571420003261,0.147619043574447];
grid on;

function power = watt2dbm(watt_value)
    % Converting power values from Watt to dbm for calculations
    %
    power = 10*log10(watt_value/(10^(-3)));
end
                  