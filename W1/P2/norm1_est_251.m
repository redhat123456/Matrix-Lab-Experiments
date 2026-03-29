function [norm_est] = norm1_est_251(B)
%% 算法 2.5.1：估计矩阵的 1 范数（优化法/盲人爬山法）
% 输入：矩阵 B
% 输出：||B||_1 的估计值

    n = size(B, 1);
    
    % 1. 选择初始向量 x，满足 ||x||_1 = 1
    x = ones(n, 1) / n;
    
    % 初始化迭代控制变量
    k_flag = 1;
    v_old = zeros(n, 1); % 用于记录上一次的符号向量以辅助判断收敛
    
    while k_flag == 1
        % 2. w = B * x
        w = B * x;
        
        % 3. v = sign(w)
        v = sign(w);
        v(v == 0) = 1; % 处理 sign(0) 的情况，约定取 1
        
        % 检查收敛：如果符号向量不再变化，说明方向已稳定
        if isequal(v, v_old)
            break;
        end
        v_old = v;
        
        % 4. z = B' * v
        z = B' * v;
        
        % 5. 检查停止准则：||z||_inf <= z' * x
        [z_inf, j] = max(abs(z));
        
        if z_inf <= (z' * x)
            k_flag = 0; % 找到局部极大值，停止迭代
        else
            % 6. 寻找下标 j 使得 |z_j| = ||z||_inf，并更新 x = e_j
            x = zeros(n, 1);
            x(j) = 1;
            k_flag = 1;
        end
    end
    
    % 7. 返回 ||w||_1
    norm_est = norm(w, 1);
end

