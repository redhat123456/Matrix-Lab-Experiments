function x = gauss_col_pivot(A, b)
    n = size(A, 1);
    for k = 1:n-1
        [~, p_rel] = max(abs(A(k:n, k)));
        p = p_rel + k - 1;
        % 同时交换 A 和 b 的行
        A([k, p], :) = A([p, k], :);
        b([k, p])    = b([p, k]);    % ← 这行原来缺失！
        
        if A(k,k) ~= 0
            A(k+1:n, k)       = A(k+1:n, k) / A(k, k);
            A(k+1:n, k+1:n)   = A(k+1:n, k+1:n) - A(k+1:n, k) * A(k, k+1:n);
            b(k+1:n)          = b(k+1:n) - A(k+1:n, k) * b(k); % ← 前代消去 b
        else
            error('矩阵奇异！');
        end
    end
    % 回代
    x = zeros(n, 1);
    x(n) = b(n) / A(n, n);
    for i = n-1:-1:1
        x(i) = (b(i) - A(i, i+1:n) * x(i+1:n)) / A(i, i);
    end
end