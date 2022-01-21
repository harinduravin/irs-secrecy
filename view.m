figure(1)
plot(L_list,mean(min_master_list));
xlabel('Number of IRS elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')


% Create plots
figure(2)
t = tiledlayout(2,2);
nexttile
plot(L_list,mean(A1_master_list));
xlabel('Number of IRS elements')
ylabel('A1 secrecy Rate(bits/sec/Hz)')
nexttile
plot(L_list,mean(A2_master_list));
xlabel('Number of IRS elements')
ylabel('A2 secrecy Rate(bits/sec/Hz)')
nexttile
plot(L_list,mean(B1_master_list));
xlabel('Number of IRS elements')
ylabel('B1 secrecy Rate(bits/sec/Hz)')
nexttile
plot(L_list,mean(B2_master_list));
xlabel('Number of IRS elements')
ylabel('B2 secrecy Rate(bits/sec/Hz)')