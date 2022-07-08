load("one_min_master_list.mat");
one_min_master_list = min_master_list;
load("two_min_master_list.mat");
two_min_master_list = min_master_list;
load("three_min_master_list.mat");
three_min_master_list = min_master_list;
load("four_min_master_list_1.mat");
four_min_master_list_1 = min_master_list; 
load("four_min_master_list_2.mat");
four_min_master_list_2 = min_master_list; 
load("five_min_master_list.mat");
five_min_master_list = min_master_list; 

L_list = 5:5:40;

plot(L_list,mean(one_min_master_list),'DisplayName','Single pair');
hold on;
plot(L_list,mean(two_min_master_list,'omitnan'),'DisplayName','Two pair');
plot(L_list,mean(three_min_master_list),'DisplayName','Three pair');
plot(L_list,mean(four_min_master_list_1),'DisplayName','Four pair');
plot(L_list,mean(five_min_master_list),'DisplayName','Five pair');
xlabel('Number of IRS elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')
legend
