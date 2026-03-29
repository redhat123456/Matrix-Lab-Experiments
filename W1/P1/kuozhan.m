%% --- 实验：随着 N 增大，分析 Hilbert 矩阵的病态性与算法表现 ---
clear; clc;

% 1. 设置参数范围
N_range = 2:2:40; % N 从 2 增加到 40
num_N = length(N_range);

% 初始化结果存储
err_results = zeros(num_N, 4); % 相对误差
res_results = zeros(num_N, 4); % 相对残差
time_results = zeros(num_N, 4); % 运行时间

% 2. 循环计算
for k = 1:num_N
    n = N_range(k);
    A = hilb(n); % 构造 Hilbert 矩阵
    y_true = ones(n, 1);
    b = A * y_true; % 构造右端项
    
    % --- 经典平方根法 ---
    tic; y1 = cholesky_solve(A, b); t1 = toc;
    % --- 算法 1.3.1 ---
    tic; y2 = cholesky_solve131(A, b); t2 = toc;
    % --- 改进平方根法 ---
    tic; y3 = Improve_cholesky_solve(A, b); t3 = toc;
    % --- MATLAB 内置 ---
    tic; y4 = A \ b; t4 = toc;
    
    % 记录时间
    time_results(k, :) = [t1, t2, t3, t4];
    
    % 记录相对误差 (norm(y-y_true)/norm(y_true))
    err_results(k, 1) = norm(y1 - y_true) / norm(y_true);
    err_results(k, 2) = norm(y2 - y_true) / norm(y_true);
    err_results(k, 3) = norm(y3 - y_true) / norm(y_true);
    err_results(k, 4) = norm(y4 - y_true) / norm(y_true);
    
    % 记录相对残差 (norm(A*y-b)/norm(b))
    res_results(k, 1) = norm(A*y1 - b) / norm(b);
    res_results(k, 2) = norm(A*y2 - b) / norm(b);
    res_results(k, 3) = norm(A*y3 - b) / norm(b);
    res_results(k, 4) = norm(A*y4 - b) / norm(b);
end

% 1. 窗口初始化
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.32]); 
methods = {'经典平方根法', '算法1.3.1', '改进平方根法(LDL^T)', 'MATLAB内置'};
colors = {'#D95319', '#EDB120', '#0072BD', '#7E2F8E'};
markers = {'o', 's', '^', 'd'};

% --- 图 (a): 相对误差 ---
ax1 = subplot(1, 3, 1); % 这里才正式创建 ax1
for i = 1:4
    semilogy(N_range, err_results(:,i), 'Marker', markers{i}, 'Color', colors{i}, 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
end
grid on; box on;
title('(a) 相对误差', 'FontName', 'SimSun', 'FontSize', 18);
xlabel('矩阵阶数 N', 'FontName', 'SimSun'); 
ylabel('误差 (对数轴)', 'FontName', 'SimSun');

% --- 图 (b): 相对残差 ---
ax2 = subplot(1, 3, 2); % 这里创建 ax2
for i = 1:4
    semilogy(N_range, res_results(:,i), 'Marker', markers{i}, 'Color', colors{i}, 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
end
grid on; box on;
title('(b) 相对残差', 'FontName', 'SimSun', 'FontSize', 18);
xlabel('矩阵阶数 N', 'FontName', 'SimSun');
ylabel('残差 (对数轴)', 'FontName', 'SimSun');

% --- 图 (c): 运行时间 ---
ax3 = subplot(1, 3, 3); % 这里创建 ax3
for i = 1:4
    plot(N_range, time_results(:,i), 'Marker', markers{i}, 'Color', colors{i}, 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
end
grid on; box on;
title('(c) 求解时间对比', 'FontName', 'SimSun', 'FontSize', 18);
xlabel('矩阵阶数 N', 'FontName', 'SimSun'); 
ylabel('时间 (s)', 'FontName', 'SimSun');

% 2. 统一美化 (必须在 subplot 之后！)
all_axes = [ax1, ax2, ax3];
for ax = all_axes
    ax.FontSize = 14;                % 刻度数字可以稍微比标题小一点，比如14
    ax.FontName = 'Times New Roman'; 
    ax.LineWidth = 1.1;           
    ax.TickDir = 'in';            
end

% 3. 图例与布局调整
lgd = legend(methods, 'Orientation', 'horizontal', 'FontName', 'SimSun', 'FontSize', 16,'Box','off');
set(lgd, 'Position', [0.1, 0.01, 0.8, 0.07], 'Units', 'normalized');
% 手动微调子图位置
set(ax1, 'Position', [0.07, 0.22, 0.25, 0.65]);
set(ax2, 'Position', [0.39, 0.22, 0.25, 0.65]);
set(ax3, 'Position', [0.71, 0.22, 0.25, 0.65]);