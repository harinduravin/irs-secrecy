load("1 pair/1_pairs_10_seeds_2_trial.mat");
one_min_master_list = min_master_list;

load("2 pairs/2_pairs_10_seeds_2_trial.mat");
two_min_master_list = min_master_list;

load("3 pairs/3_pairs_10_seeds_2_trial.mat");
three_min_master_list = min_master_list;


load("4 pair/4_pairs_10_seeds_2_trial.mat");
four_min_master_list = min_master_list; 


load("5 pairs/5_pairs_10_seeds_2_trial.mat");
five_min_master_list = min_master_list; 

load("1 pair/1_pairs_10_seeds_1_trial.mat");
one_min_master_list_e = min_master_list;

load("2 pairs/2_pairs_10_seeds_1_trial.mat");
two_min_master_list_e = min_master_list;

load("3 pairs/3_pairs_10_seeds_1_trial.mat");
three_min_master_list_e = min_master_list;


load("4 pair/4_pairs_10_seeds_1_trial.mat");
four_min_master_list_e = min_master_list; 


load("5 pairs/5_pairs_10_seeds_1_trial.mat");
five_min_master_list_e = min_master_list; 

L_list = 0:5:40;
L_list(1) = 1;

figure (1)
plot(L_list,max(mean(one_min_master_list_e),0),'DisplayName','Single pair (E1)','Marker','o','LineWidth',1.5,'Color',[0.00 0.45 0.74]);
hold on;
plot(L_list,max(mean(two_min_master_list_e),0),'DisplayName','Two pairs (E1)','Marker','*','LineWidth',1.5,'Color',[0.85 0.33 0.10]);
plot(L_list,max(mean(three_min_master_list_e),0),'DisplayName','Three pairs (E1)','Marker','^','LineWidth',1.5,'Color',[0.93 0.69 0.13]);
plot(L_list,max(mean(four_min_master_list_e),0),'DisplayName','Four pairs (E1)','Marker','x','LineWidth',1.5,'Color',[0.49 0.18 0.56]);
plot(L_list,max(mean(five_min_master_list_e),0),'DisplayName','Five pairs (E1)','Marker','+','LineWidth',1.5,'Color',[0.47 0.67 0.19]);
plot(L_list,max(mean(one_min_master_list),0),'DisplayName','Single pair (E2)','LineStyle','--','LineWidth',1.5,'Color',[0.00 0.45 0.74]);

plot(L_list,max(mean(two_min_master_list),0),'DisplayName','Two pairs (E2)','LineStyle','--','LineWidth',1.5,'Color',[0.85 0.33 0.10]);
plot(L_list,max(mean(three_min_master_list),0),'DisplayName','Three pairs (E2)','LineStyle','--','LineWidth',1.5,'Color',[0.93 0.69 0.13]);
plot(L_list,max(mean(four_min_master_list),0),'DisplayName','Four pairs (E2)','LineStyle','--','LineWidth',1.5,'Color',[0.49 0.18 0.56]);
plot(L_list,max(mean(five_min_master_list),0),'DisplayName','Five pairs (E2)','LineStyle','--','LineWidth',1.5,'Color',[0.47 0.67 0.19]);

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