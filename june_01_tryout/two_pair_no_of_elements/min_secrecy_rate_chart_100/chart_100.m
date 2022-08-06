load("one_min_master_list_1.mat");
one_min_master_list = min_master_list;
load("one_min_master_list_2.mat");
one_min_master_list = cat(1,one_min_master_list,min_master_list);

load("two_min_master_list_1.mat");
two_min_master_list = min_master_list;
% load("two_min_master_list_2.mat");
% two_min_master_list = cat(1,two_min_master_list,min_master_list);
% load("two_min_master_list_3.mat");
% two_min_master_list = cat(1,two_min_master_list,min_master_list);

load("three_min_master_list_rand.mat");
three_min_master_list_rand = min_master_list;

load("three_min_master_list_1.mat");
three_min_master_list = min_master_list;
load("three_min_master_list_2.mat");
three_min_master_list = cat(1,three_min_master_list,min_master_list);

load("four_min_master_list_1.mat");
four_min_master_list = min_master_list; 
load("four_min_master_list_2.mat");
four_min_master_list = cat(1,four_min_master_list,min_master_list);
load("four_min_master_list_3.mat");
four_min_master_list = cat(1,four_min_master_list,min_master_list);

load("five_min_master_list_1.mat");
five_min_master_list = min_master_list; 
load("five_min_master_list_2.mat");
five_min_master_list = cat(1,five_min_master_list,min_master_list);

L_list = 5:5:40;
no_irs_3_pair = zeros(1,8) + 1.1345;

plot(L_list,mean(one_min_master_list),'DisplayName','Single pair','Marker','o','LineWidth',1.5);
hold on;
plot(L_list,mean(two_min_master_list,'omitnan'),'DisplayName','Two pairs','Marker','*','LineWidth',1.5);
plot(L_list,mean(three_min_master_list),'DisplayName','Three pairs','Marker','^','LineWidth',1.5,'Color',[0.9 0.6 0]);
plot(L_list,mean(max(three_min_master_list_rand,0)),'DisplayName','Three pairs(Rnd)','LineStyle','--','LineWidth',1.5,'Color',[0.9 0.6 0]);
plot(L_list,no_irs_3_pair,'DisplayName','Three pairs(No-IRS)','LineStyle',':','LineWidth',1.5,'Color',[0.9 0.6 0]);
plot(L_list,mean(four_min_master_list),'DisplayName','Four pairs','Marker','x','LineWidth',1.5);
plot(L_list,mean(five_min_master_list),'DisplayName','Five pairs','Marker','+','LineWidth',1.5);
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

print -depsc2 chart_100_3.eps