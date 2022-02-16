fig2 = openfig('2pairs.fig');
fig3 = openfig('3pairs.fig');
fig4 = openfig('4pairs.fig');
fig5 = openfig('5pairs.fig');

axObjs2 = fig2.Children;
dataObjs2 = axObjs2.Children;
y2 = dataObjs2(1).YData;

axObjs3 = fig3.Children;
dataObjs3 = axObjs3.Children;
y3 = dataObjs3(1).YData;

axObjs4 = fig4.Children;
dataObjs4 = axObjs4.Children;
y4 = dataObjs4(1).YData;

axObjs5 = fig5.Children;
dataObjs5 = axObjs5.Children;
y5 = dataObjs5(1).YData;

close all

y_hat = [y2; y3; y4; y5];
y_matrix = y_hat';
pair_list = 2:5;

figure(1)

plot(L_list,y2,'DisplayName','2 pairs', 'Marker','diamond');
title('Minimum Secrecy Rate - Number of IRS elements')
hold on;
plot(L_list,y3,'DisplayName','3 pairs',Marker='v');
hold on;
plot(L_list,y4,'DisplayName','4 pairs',Marker='<');
hold on;
plot(L_list,y5,'DisplayName','5 pairs',Marker='>');
hold off;
lgd = legend;
lgd.FontSize = 14;
lgd.Location = 'east';
xlabel('Number of IRS elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')

figure(2)

plot(pair_list,y_matrix(3,:),'DisplayName','15 IRS elements','Marker','o');
title('Minimum secrecy rate - Number of pairs in the system')

lgd = legend;
lgd.FontSize = 14;
lgd.Location = 'east';
xlabel('Number of Pairs')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')