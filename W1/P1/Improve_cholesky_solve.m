function x = Improve_cholesky_solve(A,b)

    n = length(b);
    L = eye(n);      % 单位下三角矩阵
    d = zeros(n, 1); % 对角矩阵 D 的对角元素
 
    for j = 1:n
        % 计算 v = D(1:j-1) .* L(j,1:j-1)'
        v = d(1:j-1) .* L(j,1:j-1)';
 
        % 计算 d_j
        d(j) = A(j,j) - L(j,1:j-1) * v;
        % 计算 L(i,j), i > j
        for i = j+1:n
            L(i,j) = (A(i,j) - L(i,1:j-1) * v) / d(j);
        end
    end

    % 前代: L * y = b
    y = zeros(n, 1);
    for i = 1:n
        y(i) = b(i) - L(i,1:i-1) * y(1:i-1);
    end
 
    % 对角求解: D * z = y
    z = y ./ d;
 
    % 回代: L' * x = z
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = z(i) - L(i+1:n,i)' * x(i+1:n);
    end
end