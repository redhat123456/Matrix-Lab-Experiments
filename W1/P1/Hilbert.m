function H = Hilbert(n)
% Hilbert 生成 n x n 的 Hilbert 矩阵
% H(i,j) = 1 / (i + j - 1)
%
% 输入：
%   n - 矩阵阶数
% 输出：
%   H - n x n Hilbert 矩阵

% 初始化矩阵
H = zeros(n, n);

% 双重循环构造
for i = 1:n
    for j = 1:n
        H(i,j) = 1 / (i + j - 1);
    end
end
end