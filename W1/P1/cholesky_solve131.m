function x = cholesky_solve131(A,b)
% Cholesky 分解（平方根法，原地实现）
% 输入：A - 对称正定矩阵
% 输出：A - 下三角部分为 L，使 A = L * L'
    n = length(b);
for k = 1:n
    A(k, k) = sqrt(A(k, k)); % 计算对角元
    if k < n
        A(k+1:n, k) = A(k+1:n, k) / A(k, k); % 计算第 k 列对角线以下元素
        for j = k+1:n 
            A(j:n, j) = A(j:n, j) - A(j:n, k) * A(j, k); % 秩-1更新
        end
    end
end
    A = tril(A);
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - A(i,1:i-1) * y(1:i-1)) / A(i,i);
    end
 
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - A(i+1:n,i)' * x(i+1:n)) / A(i,i);
    end

end