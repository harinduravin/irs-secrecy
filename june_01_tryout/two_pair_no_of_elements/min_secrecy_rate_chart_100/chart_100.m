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

plot(L_list,mean(one_min_master_list),'DisplayName','Single pair','Marker','o');
hold on;
plot(L_list,mean(two_min_master_list,'omitnan'),'DisplayName','Two pair','Marker','*');
plot(L_list,mean(three_min_master_list),'DisplayName','Three pair','Marker','^');
plot(L_list,mean(four_min_master_list),'DisplayName','Four pair','Marker','x');
plot(L_list,mean(five_min_master_list),'DisplayName','Five pair','Marker','+');
ax = gca;
ax.FontSize = 12;
xlim([5 40])
xlabel('Number of IRS elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')
lgd = legend('Location','eastoutside');
lgd.FontSize = 12;
