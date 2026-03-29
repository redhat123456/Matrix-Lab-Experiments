clc; clear global; close all;
%% ========================================================
%  第二题：算法 2.5.1 估计矩阵1范数（优化法）及误差分析图表
% =========================================================

fprintf('==========================================\n');
fprintf('  第二题：1范数估计（算法2.5.1）\n');
fprintf('==========================================\n\n');

%% --- 第(1)小题：估计5到20阶Hilbert矩阵的∞范数条件数 ---
fprintf('--- 第(1)小题：Hilbert矩阵的∞范数条件数 ---\n');
fprintf('%6s %18s %18s %12s\n', '阶数n', '估计条件数(∞)', '精确条件数(∞)', '相对误差');
fprintf('%s\n', repmat('-', 1, 60));

for n = 5:20
    H = hilb(n);
    H_inv = invhilb(n);  
    
    % 注意：如果你的 norm1_est_251 被修改为返回两个值（见下文），
    % 这里只需接收第一个返回值即可。
    norm_H_inf_est   = norm1_est_251(H');
    norm_Hinv_inf_est = norm1_est_251(H_inv');

    cond_est  = norm_H_inf_est * norm_Hinv_inf_est;
    cond_exact = cond(H, inf);
    rel_err   = abs(cond_est - cond_exact) / cond_exact;

    fprintf('%6d %18.4e %18.4e %12.4e\n', n, cond_est, cond_exact, rel_err);
end
fprintf('\n');


%% --- 第(2)小题：An矩阵，估计计算解的精度 ---
fprintf('--- 第(2)小题：An矩阵 n从5到30 精度估计 ---\n');
fprintf('%6s %15s %15s %15s %15s\n', 'n', '真实误差', '估计误差上界', '真实残差', '估计/真实');
fprintf('%s\n', repmat('-', 1, 75));

% 初始化用于绘图的数组
n_vec = 5:30;
true_err_vec = zeros(length(n_vec), 1);
err_bound_vec = zeros(length(n_vec), 1);

for i = 1:length(n_vec)
    n = n_vec(i);
    An = build_An(n); % 请确保此函数存在

    rng(n);  
    x_true = randn(n, 1);
    b = An * x_true;

    x_hat = gauss_col_pivot(An, b); % 请确保此函数存在

    % 计算并保存真实相对误差
    true_err = norm(x_hat - x_true, inf) / norm(x_true, inf);
    true_err_vec(i) = true_err;

    % 估计范数
    norm_An_1   = norm1_est_251(An);
    norm_Aninv_1  = norm1_est_251(inv(An)); 
    cond_est_1 = norm_An_1 * norm_Aninv_1;

    r = b - An * x_hat;
    res = norm(r, inf) / norm(b, inf);

    % 计算并保存误差上界
    err_bound = cond_est_1 * res;
    err_bound_vec(i) = err_bound;

    fprintf('%6d %15.4e %15.4e %15.4e %15.4e\n', n, true_err, err_bound, res, err_bound/max(true_err,1e-300));
end
fprintf('\n');

%% --- 综合绘图：1范数估计性能分析 ---
% 设置阶数范围
n_range = 5:1000;
steps_record = zeros(length(n_range), 1);

% 预计算：n 增大时迭代步数的变化
for i = 1:length(n_range)
    An_temp = hilb(n_range(i));
    [~, history_temp] = norm1_est2_251(An_temp'); % 请确保函数名一致
    steps_record(i) = length(history_temp);
end

% 准备数据：n=100 时的迭代轨迹
H100 = hilb(1000);
[~, norm_history] = norm1_est2_251(H100');

% --- 布局调整 ---
% [左, 下, 宽, 高] -> 将高度从 500 调小到 400，让图看起来扁平一点
fig = figure('Name', '算法2.5.1性能分析', 'Color', 'w', 'Position', [100, 100, 1000, 500]);

% 使用 tiledlayout 更好地控制布局和全局图例
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- 左图：1范数随迭代步的变化 ---
nexttile;
p1 = plot(1:length(norm_history), norm_history, '-o', 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', '#0072BD', 'DisplayName', '1-范数估计值');
xlabel('迭代步数 (Iteration Steps)', 'FontName', 'SimSun', 'FontSize', 16);
ylabel('1-范数估计值', 'FontName', 'SimSun', 'FontSize', 16);
title('1-范数随迭代步的变化 (n=1000)', 'FontName', 'SimSun', 'FontSize', 16, 'FontWeight', 'bold');
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% --- 右图：迭代步数随 n 的变化 ---
nexttile;
p2 = plot(n_range, steps_record, '-s', 'LineWidth', 1.5, ...
    'Color', '#D95319', 'MarkerFaceColor', '#D95319', 'DisplayName', '所需迭代步数');
xlabel('矩阵阶数 (Matrix Order n)', 'FontName', 'SimSun', 'FontSize', 16);
ylabel('迭代次数 (Iteration Count)', 'FontName', 'SimSun', 'FontSize', 16);
title('迭代步数随矩阵阶数 n 的变化', 'FontName', 'SimSun', 'FontSize', 16, 'FontWeight', 'bold');
ylim([0, max(steps_record)+2]); 
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% --- 全局图例设置 ---
% 将图例绑定到布局 t 上，而不是单个子图
lgd = legend([p1, p2], 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % 关键：置于布局的最下方
set(lgd, 'FontName', 'SimSun', 'FontSize', 16);



%% --- 综合绘图：算法2.5.1 对 An 矩阵的性能与误差全景分析 ---

% 1. 数据准备
n_range = 5:40; 
steps_record = zeros(length(n_range), 1);
cond_est_vec = zeros(length(n_range), 1);
cond_real_vec = zeros(length(n_range), 1);
true_err_vec = zeros(length(n_range), 1);  % 存储真实相对误差
err_bound_vec = zeros(length(n_range), 1); % 存储估计误差上界

for i = 1:length(n_range)
    n = n_range(i);
    An = build_An(n); % 调用你的 build_An 函数
    
    % --- 求解线性方程组 (用于计算真实误差) ---
    rng(n); 
    x_true = randn(n, 1);
    b = An * x_true;
    x_hat = gauss_col_pivot(An, b); % 调用你的高斯消去函数
    
    % 真实误差 (inf范数)
    true_err_vec(i) = norm(x_hat - x_true, inf) / norm(x_true, inf);
    
    % --- 算法 2.5.1 估计 1-范数条件数 ---
    [norm_An, hist] = norm1_est2_251(An);
    [norm_An_inv, ~] = norm1_est2_251(inv(An));
    
    cond_est_vec(i) = norm_An * norm_An_inv;
    cond_real_vec(i) = cond(An, 1);
    steps_record(i) = length(hist);
    
    % 估计误差上界: cond(A) * (||r||/||b||)
    r = b - An * x_hat;
    err_bound_vec(i) = cond_est_vec(i) * (norm(r, inf) / norm(b, inf));
end

% 准备左图数据 (以最大规模 n=30 为例观察单次迭代)
An_max = build_An(max(n_range));
[~, norm_history] = norm1_est2_251(An_max);

% 2. 创建大图 (高度 380 使其扁平)
fig = figure('Name', 'An矩阵数值实验分析', 'Color', 'w', 'Position', [50, 200, 1300, 450]);
t = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

% --- 图1：1-范数随迭代步的变化 (算法收敛性) ---
nexttile;
p1 = plot(1:length(norm_history), norm_history, '-o', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'MarkerFaceColor', '#0072BD', 'DisplayName', '1-范数估计值轨迹');
xlabel('迭代步数 (Steps)', 'FontName', 'SimSun', 'FontSize', 16);
ylabel('估计值', 'FontName', 'SimSun', 'FontSize', 16);
title('1-范数迭代轨迹 (n=50)', 'FontName', 'SimSun', 'FontSize', 15, 'FontWeight', 'bold');
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% --- 图2：条件数对比 (验证估计精度) ---
nexttile;
p2 = semilogy(n_range, cond_real_vec, 'k-', 'LineWidth', 2, 'DisplayName', '理论条件数 \kappa_1(A_n)');
hold on;
p3 = semilogy(n_range, cond_est_vec, 'r--', 'LineWidth', 1.5, 'DisplayName', '估计条件数 \kappa_{est}');
xlabel('矩阵阶数 n', 'FontName', 'SimSun', 'FontSize', 16);
ylabel('条件数', 'FontName', 'SimSun', 'FontSize', 16);
title('条件数估计精度对比', 'FontName', 'SimSun', 'FontSize', 16, 'FontWeight', 'bold');
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% --- 图3：相对误差与误差上界 (精度分析) ---
nexttile;
p4 = semilogy(n_range, true_err_vec, 'b-d', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '实际相对误差');
hold on;
p5 = semilogy(n_range, err_bound_vec, 'm-s', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '估计误差上界');
xlabel('矩阵阶数 n', 'FontName', 'SimSun', 'FontSize', 16);
ylabel('误差水平', 'FontName', 'SimSun', 'FontSize', 16);
title('解的精度与误差评估', 'FontName', 'SimSun', 'FontSize', 16, 'FontWeight', 'bold');
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% --- 全局图例设置：放置于正下方 ---
lgd = legend([p1, p2, p3, p4, p5], 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; 
set(lgd, 'FontName', 'SimSun', 'FontSize', 16, 'Box', 'on', 'EdgeColor', [.5 .5 .5]);

% 保持美观的缩放
set(gcf, 'Units', 'normalized');