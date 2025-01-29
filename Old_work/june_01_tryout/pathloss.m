f = 1;
d = 1:5:201;
pl1 = -28-20*log10(f)-22*log10(d);
pl2 = -22.7-26*log10(f)-36.7*log10(d);

% pl3 = -32.8-20*log10(f)-16.9*log10(d);
% pl4 = -11.5-20*log10(f)-43.3*log10(d);

pl3 = -15-20*log10(f)-20*log10(d);
pl4 = -15-20*log10(f)-30*log10(d);

figure(1)
plot(d,pl1,'-<',d,pl2,'->',d,pl3,'-^',d,pl4,'-o');
% plot(d,pl3,'-<',d,pl4,'->');
% xlabel('Number of elements')
% ylabel('Minimum secrecy Rate(bits/sec/Hz)')