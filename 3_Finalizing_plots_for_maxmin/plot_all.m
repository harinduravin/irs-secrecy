clear;
% close all;
load('config2.mat');

maxmin_array = zeros(1,length(L_set));
rand_data_min_array = zeros(1,length(L_set)); 
rand_data_max_array = zeros(1,length(L_set)); 
rand_data_popt_array = zeros(1,length(L_set)); 

p1_array = zeros(1,length(L_set));
p2_array = zeros(1,length(L_set));
p1_a_array = zeros(1,length(L_set));
p1_b_array = zeros(1,length(L_set));
p1_c_array = zeros(1,length(L_set));
p2_a_array = zeros(1,length(L_set));
p2_b_array = zeros(1,length(L_set));
p2_c_array = zeros(1,length(L_set));

for i = 1:length(L_set)
    maxmin_array(i) = mean(store_data(i,:,1),"omitnan");
    rand_data_min_array(i) = mean(rand_data_min(i,:)); 
    rand_data_max_array(i) = mean(rand_data_max(i,:));
    rand_data_popt_array(i) = mean(rand_data_popt(i,:));
    p1_array(i) = mean(store_data(i,:,2));
    p2_array(i) = mean(store_data(i,:,3));
    p1_a_array(i) = mean(store_data(i,:,4));
    p1_b_array(i) = mean(store_data(i,:,5));
    p1_c_array(i) = mean(store_data(i,:,6));
    p2_a_array(i) = mean(store_data(i,:,7));
    p2_b_array(i) = mean(store_data(i,:,8));
    p2_c_array(i) = mean(store_data(i,:,9));
end

no_irs_popt_array = mean(store_data(9,:,1));
no_irs_min_array = mean(rand_data_min(9,:));
no_irs_max_array = mean(rand_data_max(9,:));

figure
plot(L_set,maxmin_array,'-o','DisplayName','Power and IRS optimized'); hold on;
plot(L_set,rand_data_min_array,'-v','DisplayName','Minimum power for all with random IRS');
plot(L_set,rand_data_max_array,'-<','DisplayName','Maximum power for all with random IRS');
plot(L_set,rand_data_popt_array,'-o','DisplayName','Power optimized with random IRS');
yline(no_irs_popt_array,'DisplayName','No-IRS Power Opt','LineStyle','-.')
yline(no_irs_min_array,'DisplayName','No-IRS Power min')
yline(no_irs_max_array,'DisplayName','No-IRS Power max')
grid on;
legend('Location','southwest')
xlabel('Number of IRS elements')
ylabel('Max-min Secrecy Rate (bps/Hz)')
% ylim([0 2.5])
set(gca,'GridLineStyle','--')

figure
plot(L_set,maxmin_array,'-o','DisplayName','maxmin'); hold on;
plot(L_set,p1_array,'-v','DisplayName','p1');
plot(L_set,p2_array,'-<','DisplayName','p2');
plot(L_set,p1_a_array,'-o','DisplayName','p1a');
plot(L_set,p1_b_array,'-o','DisplayName','p1b');
plot(L_set,p1_c_array,'-o','DisplayName','p1c');
plot(L_set,p2_a_array,'-o','DisplayName','p2a');
plot(L_set,p2_b_array,'-o','DisplayName','p2b');
plot(L_set,p2_c_array,'-o','DisplayName','p2c');
% yline(no_irs_popt_array,'DisplayName','No-IRS Power Opt','LineStyle','-.')
% yline(no_irs_min_array,'DisplayName','No-IRS Power min')
% yline(no_irs_max_array,'DisplayName','No-IRS Power max')
grid on;
legend('Location','southwest')
xlabel('Number of IRS elements')
ylabel('Max-min Secrecy Rate (bps/Hz)')
% ylim([0 2.5])
set(gca,'GridLineStyle','--')
