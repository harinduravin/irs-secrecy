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

figure(1)
plot(L_list,y2);
hold on;
plot(L_list,y3);
hold on;
plot(L_list,y4);
hold on;
plot(L_list,y5);
hold on;
xlabel('Number of elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')