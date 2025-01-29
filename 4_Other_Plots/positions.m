coords = [-20 -10; -20 -40; 20 -10; 20 -40; 0 10; 0 5;0 0; -15 -35; -15 -45];

coords_shape = size(coords);
n = (coords_shape(1) - 5) / 2;

a = coords(:, 1);
b = coords(:, 2);
c = {};
for i = 1:n
    c{end + 1} = ['$A_', num2str(i), '$'];
    c{end + 1} = ['$B_', num2str(i), '$'];
end

c{end + 1} = '$V$';
c{end + 1} = '$W$';
c{end + 1} = '$X$';
c{end + 1} = '$Y$';
c{end + 1} = '$Z$';

sz = 70;
leg_a = a(1:2 * n);
leg_b = b(1:2 * n);

% Label will be "North" of the datapoint with 0.1 spacing
figure(3)
hold on;

% Adding scatter plots and dynamic legend labels
for i = 1:n
    scatter(leg_a(2*i), leg_b(2*i), sz, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'DisplayName', ['$B_', num2str(i), '$ (', num2str(leg_a(2*i)), ', ', num2str(leg_b(2*i)), ')'])
    scatter(leg_a(2*i-1), leg_b(2*i-1), sz, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'DisplayName', ['$A_', num2str(i), '$ (', num2str(leg_a(2*i-1)), ', ', num2str(leg_b(2*i-1)), ')'])
end

scatter(a(2 * n + 1), b(2 * n + 1), sz, 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'DisplayName', ['$V$ (', num2str(a(2 * n + 1)), ', ', num2str(b(2 * n + 1)), ')'])
scatter(a(2 * n + 2), b(2 * n + 2), sz, 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'DisplayName', ['$W$ (', num2str(a(2 * n + 2)), ', ', num2str(b(2 * n + 2)), ')'])
scatter(a(2 * n + 3), b(2 * n + 3), sz, 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'DisplayName', ['$X$ (', num2str(a(2 * n + 3)), ', ', num2str(b(2 * n + 3)), ')'])
scatter(a(2 * n + 4), b(2 * n + 4), sz, 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'DisplayName', ['$Y$ (', num2str(a(2 * n + 4)), ', ', num2str(b(2 * n + 4)), ')'])
scatter(a(2 * n + 5), b(2 * n + 5), sz, 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'DisplayName', ['$Z$ (', num2str(a(2 * n + 5)), ', ', num2str(b(2 * n + 5)), ')'])

% Adding labels close to each point
for i = 1:length(a)
    text(a(i) + 1.3, b(i) + 1.3, c{i}, 'FontSize', 12, 'FontName', 'Times', 'Interpreter', 'latex')
end

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
ax.Box = 'on';

xlabel('X (m)', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12)
ylabel('Y (m)', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12)

lgd = legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex', 'FontName', 'Times');
grid on;
axis equal;
ylim([-50 15])
xlim([-45 100])
