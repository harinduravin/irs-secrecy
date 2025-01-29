clear;

load('convergence_sum_power_3.mat','fractional_increases','diffwithexh','num_iterations','obj_evolution','optimum','numrealizations')

% Identify columns that contain zeros
columns_with_zeros = any(obj_evolution <= 0, 1);

% Remove columns that contain zeros
obj_evolution(:, columns_with_zeros) = [];
optimum(columns_with_zeros) = [];

fractional_increases = zeros(num_iterations-1,numrealizations);
diffwithexh = zeros(num_iterations,numrealizations);
for r = 1:size(obj_evolution,2)
    fractional_increases(:,r) = diff(obj_evolution(:,r))./obj_evolution(1:end-1,r);
    diffwithexh(:,r) = abs(obj_evolution(:,r) - optimum(r))/optimum(r);
end

fractional_increases_sum = fractional_increases;
diffwithexh_sum = diffwithexh;
num_iterations_sum = num_iterations;

load('converge_power_2.mat')

% Identify columns that contain zeros
columns_with_zeros = any(obj_evolutionset == 0, 1);

% Remove columns that contain zeros
obj_evolutionset(:, columns_with_zeros) = [];
optvalueexh(columns_with_zeros) = [];

fractional_increases = zeros(num_iterations-1,numrealizations);
diffwithexh = zeros(num_iterations,numrealizations);
for r = 1:size(obj_evolutionset,2)
    fractional_increases(:,r) = diff(obj_evolutionset(:,r))./obj_evolutionset(1:end-1,r);
    diffwithexh(:,r) = abs(obj_evolutionset(:,r) - optvalueexh(r))/optvalueexh(r);
end

% Set the colors
color1 = [1, 0, 0]; % Red
color2 = [0, 0, 1]; % Blue
color3 = [0, 0.5, 0]; % Dark Green
color4 = [0.5, 0, 0.5]; % Dark Purple
color5 = [0.5, 0.25, 0]; % Dark Orange


% Set the marker size
marker_size = 4; % Adjust the marker size as needed

figure;

subplot(1,2,1)

semilogy(1:num_iterations, mean(diffwithexh,2),'-v','Color',color1, 'MarkerSize', marker_size, 'MarkerFaceColor', color1, DisplayName='$\Delta_i^{Alt}$ (N = 2)')
hold on;
semilogy(1:num_iterations-1, mean(fractional_increases,2),'-o','Color',color1, 'MarkerSize', marker_size, 'MarkerFaceColor', color1, DisplayName='$\Delta_i$ (N = 2)')
semilogy(1:num_iterations_sum, mean(diffwithexh_sum,2),'-v','Color',color2, 'MarkerSize', marker_size, 'MarkerFaceColor', color2, DisplayName='$\Delta_i^{Alt}$ (N = 1)');
semilogy(1:num_iterations_sum-1, mean(abs(fractional_increases_sum),2),'-o','Color',color2, 'MarkerSize', marker_size, 'MarkerFaceColor', color2, DisplayName='$\Delta_i$ (N = 1)')

ylim([10^(-6) 10^(1)])
xlim([0 40])
xticks([0 10 20 30 40])
yticks(10.^[-6:1])
xlabel({'Number of Iterations','(a)'},'FontSize',12,Interpreter='latex')
ylabel('Fractional Increase','FontSize',12,Interpreter='latex')
grid on;
legend('FontSize',12,Location="northeast",Interpreter='latex');
set(gca,'FontSize',12)
ax = gca;
ax.YMinorGrid = 'off';
ax.GridLineStyle = '--';

load('power_conv_variables.mat','optPa','optPb2','optPb','optPa2','optvalueexh','obj_evolutionset','num_iterations','Pa_evolutionset','Pb_evolutionset','Pa2_evolutionset','Pb2_evolutionset','r')

subplot(1,2,2)
yline(optPa(r)/optPb2(r), '--', 'LabelHorizontalAlignment', 'left','Color',color3,'HandleVisibility','off');
hold on;
yline(optPb(r)/optPb2(r), '--','HandleVisibility','off', 'LabelHorizontalAlignment', 'left','Color',color4);
yline(optPa2(r)/optPb2(r), '--','HandleVisibility','off', 'LabelHorizontalAlignment', 'left','Color',color5);
yline(optvalueexh(r), '--','HandleVisibility','off', 'LabelHorizontalAlignment', 'left','Color',[0,0,0]);
% yline(optPb2/optPb2, 'g--', 'Pb2opt', 'LabelHorizontalAlignment', 'left');
plot(1:num_iterations, Pa_evolutionset(:,r)./Pb2_evolutionset(:,r), '-o','DisplayName','$P_{A1}/P_{B2}$','Color',color3, 'MarkerSize', marker_size, 'MarkerFaceColor', color3);
plot(1:num_iterations, Pb_evolutionset(:,r)./Pb2_evolutionset(:,r), '-o','DisplayName','$P_{B1}/P_{B2}$','Color',color4, 'MarkerSize', marker_size, 'MarkerFaceColor', color4);
plot(1:num_iterations, Pa2_evolutionset(:,r)./Pb2_evolutionset(:,r), '-o','DisplayName','$P_{A2}/P_{B2}$','Color',color5, 'MarkerSize', marker_size, 'MarkerFaceColor', color5);
plot(1:num_iterations, obj_evolutionset(:,r), '-o','DisplayName','Objective','Color',[0,0,0], 'MarkerSize', marker_size, 'MarkerFaceColor', [0,0,0]);
% plot(1:num_iterations, Pb2_evolution, 'b-o');
% ylim([Pmin, Pmax]);
ylim([0.5, 2.3]);
xlim([0, num_iterations+1]);
xlabel({'Number of Iterations','(b)'},Interpreter='latex');
ylabel('Objective',Interpreter='latex');
xticks([0 10 20 30 40])
% title('Optimization Results: Paopt and Pbopt');
% legend('Paopt', 'Pbopt', 'Location', 'best');
legend('FontSize',12,'location','northeast','Interpreter','latex');
grid on;
box on;
set(gca,'FontSize',12)
ax = gca;
ax.GridLineStyle = '--';
hold off;