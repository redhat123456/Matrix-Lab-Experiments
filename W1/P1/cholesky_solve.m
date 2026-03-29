function x = cholesky_solve(A, b)
%% 函数1：平方根法（Cholesky分解）求解 Ax = b
    n = length(b);
    L = zeros(n, n);
     for j = 1:n
        s = A(j,j) - sum(L(j,1:j-1).^2);
        L(j,j) = sqrt(s);
 % 计算下三角元素
        for i = j+1:n
            L(i,j) = (A(i,j) - sum(L(i,1:j-1).*L(j,1:j-1))) / L(j,j);
        end
    end
 
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - L(i,1:i-1) * y(1:i-1)) / L(i,i);
    end
 
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - L(i+1:n,i)' * x(i+1:n)) / L(i,i);
    end
end