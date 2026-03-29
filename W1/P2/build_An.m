function An = build_An(n)
    % 1. 构造单位阵 I (对角线为 1)
    I = eye(n);
    
    % 2. 构造严格下三角阵 L (元素全为 -1)
    L = -tril(ones(n), -1);
    
    % 3. 构造末列修正阵 
    C = zeros(n);
    C(1:n-1, n) = 1;

    An = I + L + C;
end