load("1 pair/1_pairs_10_seeds_1_trial.mat");
one_min_master_list = min_master_list;

load("2 pairs/2_pairs_10_seeds_1_trial.mat");
two_min_master_list = min_master_list;

load("3 pairs/3_pairs_10_seeds_1_trial.mat");
three_min_master_list = min_master_list;


load("4 pair/4_pairs_10_seeds_1_trial.mat");
four_min_master_list = min_master_list; 


load("5 pairs/5_pairs_10_seeds_1_trial.mat");
five_min_master_list = min_master_list; 

L_list = 0:5:40;
L_list(1) = 1;

figure (1)
plot(L_list,max(mean(one_min_master_list),0),'DisplayName','Single pair','Marker','o','LineWidth',1.5);
hold on;
plot(L_list,max(mean(two_min_master_list),0),'DisplayName','Two pairs','Marker','*','LineWidth',1.5);
plot(L_list,max(mean(three_min_master_list),0),'DisplayName','Three pairs','Marker','^','LineWidth',1.5,'Color',[0.9 0.6 0]);
% plot(L_list,mean(max(three_min_master_list_rand,0)),'DisplayName','Three pairs(Rnd)','LineStyle','--','LineWidth',1.5,'Color',[0.9 0.6 0]);
% plot(L_list,no_irs_3_pair,'DisplayName','Three pairs(No-IRS)','LineStyle',':','LineWidth',1.5,'Color',[0.9 0.6 0]);
plot(L_list,max(mean(four_min_master_list),0),'DisplayName','Four pairs','Marker','x','LineWidth',1.5);
plot(L_list,max(mean(five_min_master_list),0),'DisplayName','Five pairs','Marker','+','LineWidth',1.5);
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
xlim([5 40])
xlabel('Number of IRS elements','FontName','Times','FontSize',12)
ylabel('Minimum secrecy Rate(bits/sec/Hz)','FontName','Times','FontSize',12)
lgd = legend('Location','best');
lgd.FontSize = 10;
lgd.FontName = 'Times';
lgd.Position = [0.6656    0.3651    0.2245    0.2429];
grid on;
title({'Variation of Minimum secrecy rate with',' Number of IRS elements'})

% print -depsc2 chart_100_3.eps