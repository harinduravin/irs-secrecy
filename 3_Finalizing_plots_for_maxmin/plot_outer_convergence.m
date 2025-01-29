load("config4.mat")

altcolors = ['#1f77b4'; '#ff7f0e'; '#2ca02c'; '#d62728'; '#9467bd'; '#8c564b'; '#e377c2'; '#7f7f7f'; ...
          '#bcbd22'; '#17becf'; '#aec7e8'; '#ffbb78'; '#98df8a'; '#ff9896'; '#c5b0d5'; '#c49c94'; ...
          '#f7b6d2'; '#c7c7c7'; '#dbdb8d'; '#9edae5'; '#ff9896'; '#c5b0d5'; '#c49c94'; '#f7b6d2'; ...
          '#c7c7c7'; '#dbdb8d'; '#9edae5'; '#17becf'; '#aec7e8'; '#ffbb78'; '#98df8a'; '#ff7f0e'];


% r = 30;
num_iterations = 4;
numrealizations = R;
fractional_increases = zeros(8,num_iterations-1,numrealizations);

for l = 1:8
    % diffwithexh = zeros(num_iterations,numrealizations);
    obj_evolution = zeros(num_iterations,numrealizations);
    for r = 1:size(convergence_data,2)
        obj_evolution(:,r) = [rand_data_min(l,r) squeeze(convergence_data(l,r,:))'];
    end
    
    % Identify columns that contain zeros
    columns_with_zeros = any(obj_evolution == 0, 1);
    
    % Remove columns that contain zeros
    obj_evolution(:, columns_with_zeros) = [];

    for r = 1:size(obj_evolution,2)
        fractional_increases(l,:,r) = diff(obj_evolution(:,r))./obj_evolution(1:end-1,r);
    end
end

% Set the colors
color1 = [1, 0, 0]; % Red
color2 = [0, 0, 1]; % Blue

% Set the marker size
marker_size = 4; % Adjust the marker size as needed

figure;

% semilogy(1:num_iterations, mean(diffwithexh,2),'-v','Color',color1, 'MarkerSize', marker_size, 'MarkerFaceColor', color1, DisplayName='N = 2, Exh')
% hold on;

for l = 1:8
    semilogy(1:num_iterations-1, mean(squeeze(fractional_increases(l,:,:)),2),'-o', 'Color',altcolors(l,:), 'MarkerSize', marker_size, 'MarkerFaceColor', altcolors(l,:),  DisplayName=['L = ', num2str(l*5)]);hold on;
end

ylim([10^(-4) 10])
xlim([1 3])
xticks([1 2 3])
xlabel('Number of Iterations','FontSize',14,Interpreter='latex')
ylabel('Fractional Increase','FontSize',14,Interpreter='latex')
grid on;
legend('FontSize',12','NumColumns',2,Location="southwest",Interpreter='latex');
set(gca,'FontSize',12)
ax = gca;
ax.YMinorGrid = 'off';
ax.GridLineStyle = '--';