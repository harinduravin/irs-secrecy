% Produce the data
% coords = [87.49 82.30;64.03 64.11;-45.12 65.73;
% -58.74 84.97;-86.17 3.14;-64.26 -13.22;
% 10.86 20.34;-15.31 9.52;61.38 -5.16;
% 91.09 12.03;0.00 97.50;0.00 -97.50;
% ];

% coords = [50.13 76.03;31.56 40.60;-85.14 -38.49;
% -88.72 -78.33;-51.19 61.90;-11.20 61.29;
% 0.00 97.50;0.00 -97.50];

% coords = [50.13 76.03;31.56 40.60;-85.14 -38.49;
% -88.72 -78.33;-51.19 61.90;-11.20 61.29;
% 77.36 35.80;89.39 73.95;-60.15 7.59;
% -23.08 -7.43;0.00 97.50;0.00 -97.50;
% ];

% coords = [50.13 76.03;31.56 40.60;-85.14 -38.49;
% -88.72 -78.33;-51.19 61.90;-11.20 61.29;
% 77.36 35.80;89.39 73.95;-60.15 7.59;
% -23.08 -7.43;0.00 97.50;0.00 -97.50;
% ];

% coords = [-15.67 -49.65;-51.77 -32.41;41.93 75.78;
% 34.43 36.49;-70.43 71.51;-32.27 83.52;
% 0.00 97.50;100.00 0.00;
% ];

% coords = [50.13 76.03;31.56 40.60;-85.14 -38.49;
% -88.72 -78.33;-51.19 61.90;-11.20 61.29;
% 77.36 35.80;89.39 73.95;-60.15 7.59;
% -23.08 -7.43;0.00 97.50;0.00 -60;
% ];

coords = [50.13 76.03;31.56 40.60;-85.14 -38.49;
-88.72 -78.33;-51.19 61.90;-11.20 61.29;
0.00 97.50;0.00 0.00;
];

coords = [50.13 76.03;31.56 40.60;-85.14 -38.49;
-88.72 -78.33;
0.00 97.50;0.00 0.00];

coords_shape = size(coords);
n = (coords_shape(1) -2)/2;

a=coords(:,1);
b=coords(:,2);
c = {};
for i = 1:n
    c{end+1} = 'A'+string(i);
    c{end+1} = 'B'+string(i);
end
c{end+1} = 'IRS';
c{end+1} = 'E';

sz = 70;
leg_a = a(1:2*n);
leg_b = b(1:2*n);

% Label will be "North" of the datapoint with 0.1 spacing
figure(3)
scatter(leg_a(2:2:end),leg_b(2:2:end),sz,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1],'DisplayName','Group A')
hold on;
scatter(leg_a(1:2:end),leg_b(1:2:end),sz,'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1],'DisplayName','Group B')
scatter(a(2*n+1),b(2*n+1),sz,'s','MarkerFaceColor',[0.3 0.5 0.8],'MarkerEdgeColor',[1 1 1],'DisplayName','IRS')
scatter(a(2*n+2),b(2*n+2),sz,'p','MarkerFaceColor',[0.8 0.5 0.3],'MarkerEdgeColor',[0 0.5 0],'DisplayName' ,'Eavesdropper')
labelpoints(a,b,c,'N',0.2,'FontSize',10,'FontName','Times')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
ax.Box = 'on';

ylim([-100 100])
xlim([-100 100])
xlabel('X (m)','Interpreter','latex','FontName','Times','FontSize',12)
ylabel('Y (m)', 'Interpreter','latex','FontName','Times','FontSize',12)
xticks([-100 -75 -50 -25 0 25 50 75 100])
yticks([-100 -75 -50 -25 0 25 50 75 100])
% lgd = legend('Location','eastoutside');
% lgd.FontSize = 12;
lgd = legend('Location','southeast','FontSize',12,'Interpreter','latex',FontName = 'Times');
grid on;
axis equal;

print -depsc2 positions.eps