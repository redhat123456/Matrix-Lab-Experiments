function [est_norm, history] = norm1_est2_251(A)
% norm1_est_251: 估计矩阵 A 的 1-范数 (基于优化法/Hager-Higham算法)
% 输入:
%   A       - 待估计的矩阵
% 输出:
%   est_norm- 估计得到的 1-范数近似值
%   history - 记录每次迭代中 1-范数估计值的列向量，用于绘制收敛曲线

    n = size(A, 2);
    
    % 1. 初始化
    % 初始测试向量 x 通常取各分量均为 1/n 的向量
    x = ones(n, 1) / n; 
    
    history = [];    % 初始化历史记录数组
    max_iter = 10;   % 设置最大迭代次数，该算法通常在 4~5 步内收敛
    
    for iter = 1:max_iter
        % 计算 y = Ax
        y = A * x;
        
        % 当前的 1-范数估计值
        current_est = norm(y, 1);
        
        % ★ 核心修改：将当前的估计值存入历史记录
        history = [history; current_est];
        
        % 计算 z = A^T * sign(y)
        % 注意：为了防止 y 中有 0 导致 sign 变为 0，通常将 0 替换为 1
        s = sign(y);
        s(s == 0) = 1; 
        z = A' * s;
        
        % 找到 z 中绝对值最大的元素的索引 j
        [max_z, j] = max(abs(z));
        
        % 2. 停止准则
        % 如果 z 的最大无穷范数 <= z^T * x，说明已经达到局部最大值
        if max_z <= z' * x
            break;
        end
        
        % 3. 更新 x
        % 构造新的单位向量 e_j
        x_new = zeros(n, 1);
        x_new(j) = 1;
        
        % 如果新选出的单位向量和上一步完全一样，也会陷入死循环，应提前跳出
        if isequal(x, x_new)
            break;
        end
        
        % 将 x 更新为 e_j，进入下一次迭代
        x = x_new;
    end
    
    % 最终的估计结果就是历史记录里的最后一个值
    est_norm = history(end);
end