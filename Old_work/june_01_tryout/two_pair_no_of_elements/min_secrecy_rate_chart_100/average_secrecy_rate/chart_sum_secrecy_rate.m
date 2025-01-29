load("one_all_master_list_1.mat");
one_all_master_list = all_master_list;
load("one_all_master_list_2.mat");
one_all_master_list = cat(3,one_all_master_list,all_master_list);
sum_one_all_master_list = sum(mean(one_all_master_list,3),2);

load("two_all_master_list_1.mat");
two_all_master_list = all_master_list;
load("two_all_master_list_2.mat");
two_all_master_list = cat(3,two_all_master_list,all_master_list);
load("two_all_master_list_3.mat");
two_all_master_list = cat(3,two_all_master_list,all_master_list);
sum_two_all_master_list = sum(mean(two_all_master_list,3),2);

load("three_all_master_list_1.mat");
three_all_master_list = all_master_list;
load("three_all_master_list_2.mat");
three_all_master_list = cat(3,three_all_master_list,all_master_list);
sum_three_all_master_list = sum(mean(three_all_master_list,3),2);

load("four_all_master_list_1.mat");
four_all_master_list = all_master_list;
load("four_all_master_list_2.mat");
four_all_master_list = cat(3,four_all_master_list,all_master_list);
load("four_all_master_list_3.mat");
four_all_master_list = cat(3,four_all_master_list,all_master_list);
sum_four_all_master_list = sum(mean(four_all_master_list,3),2);

load("five_all_master_list_1.mat");
five_all_master_list = all_master_list;
load("five_all_master_list_2.mat");
five_all_master_list = cat(3,five_all_master_list,all_master_list);
sum_five_all_master_list = sum(mean(five_all_master_list,3),2);


all_data = [sum_one_all_master_list';sum_two_all_master_list';sum_three_all_master_list';sum_four_all_master_list'; sum_five_all_master_list']';
L_list = 1:5;

plot(L_list,all_data(1,:),'DisplayName','5 IRS elements','Marker','o');
hold on;
plot(L_list,all_data(2,:),'DisplayName','10 IRS elements','Marker','*');
plot(L_list,all_data(3,:),'DisplayName','15 IRS elements','Marker','^');
plot(L_list,all_data(4,:),'DisplayName','20 IRS elements','Marker','x');
plot(L_list,all_data(5,:),'DisplayName','25 IRS elements','Marker','+');
plot(L_list,all_data(6,:),'DisplayName','30 IRS elements','Marker','<');
plot(L_list,all_data(7,:),'DisplayName','35 IRS elements','Marker','>');
plot(L_list,all_data(8,:),'DisplayName','40 IRS elements','Marker','s');
ax = gca;
ax.FontSize = 12;
xlabel('Number of Pairs')
ylabel('Sum of individual secrecy rates(bits/sec/Hz)')
lgd = legend('Location','eastoutside');
lgd.FontSize = 12;
