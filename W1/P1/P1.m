clc; clear global; close all;

fprintf('==========================================\n');
fprintf('  第一题: 平方根法与改进的平方根法\n');
fprintf('==========================================\n\n');

fprintf('--- 第(1)小题: 100阶矩阵 ---\n');
n1 = 100;

A1 = 10 * eye(n1) + diag(ones(n1-1,1), 1) + diag(ones(n1-1,1), -1);

rng(42);
b1 = randn(n1, 1);

tic;
x1_chol = cholesky_solve(A1, b1);
time_chol = toc;

tic;
x2_chol = cholesky_solve131(A1, b1);
time_ldlt = toc;

tic;
x1_ldlt = Improve_cholesky_solve(A1, b1);
time_ildlt = toc;

tic;
x1_ref = A1 \ b1;
time_mat = toc;

fprintf('\n====== 第(1)小题 求解结果 ======\n');
disp('经典平方根法求解结果 x:');
disp(x1_chol');
disp('算法1.3.1的平方根法求解结果 x:');
disp(x2_chol');
disp('改进平方根法求解结果 x:');
disp(x1_ldlt');
disp('MATLAB内置法求解结果 x:');
disp(x1_ref');

err1_chol = norm(x1_chol - x1_ref) / norm(x1_ref);
err2_chol = norm(x2_chol - x1_ref) / norm(x1_ref);
err1_ldlt = norm(x1_ldlt - x1_ref) / norm(x1_ref);

res1_chol = norm(A1 * x1_chol - b1) / norm(b1);
res2_chol = norm(A1 * x2_chol - b1) / norm(b1);
res1_ldlt = norm(A1 * x1_ldlt - b1) / norm(b1);

fprintf('\n====== 相对误差对比 ======\n');
fprintf('经典的平方根法  相对误差: %.4e,  残差: %.4e\n', err1_chol, res1_chol);
fprintf('算法1.3.1的平方根法        相对误差: %.4e,  残差: %.4e\n', err2_chol, res2_chol);
fprintf('改进平方根法    相对误差: %.4e,  残差: %.4e\n\n', err1_ldlt, res1_ldlt);

fprintf('\n====== 求解时间对比 ======\n');
fprintf('经典的平方根法: %10.6f 秒\n', time_chol);
fprintf('算法1.3.1的平方根法: %10.6f 秒\n', time_ldlt);
fprintf('改进平方根法: %10.6f 秒\n', time_ildlt);
fprintf('MATLAB内置法: %10.6f 秒\n', time_mat);
fprintf('==========================\n');

fprintf('--- 第(2)小题: 40阶 Hilbert 矩阵 ---\n');
n2 = 40;

A2 = Hilbert(n2);

b2 = zeros(n2, 1);
for i = 1:n2
    b2(i) = sum(1 ./ (i + (1:n2) - 1));
end

y_true = ones(n2, 1);

tic;
y1_chol = cholesky_solve(A2, b2);
time_chol2 = toc;

tic;
y2_chol = cholesky_solve131(A2, b2);
time_ldlt2 = toc;

tic;
y2_ldlt = Improve_cholesky_solve(A2, b2);
time_ildlt2 = toc;

tic;
y2_ref = A2 \ b2;
time_mat2 = toc;

fprintf('\n====== 第(2)小题 求解结果 ======\n');
disp('经典平方根法求解结果 y:');
disp(y1_chol');
disp('算法1.3.1的平方根法求解结果 y:');
disp(y2_chol');
disp('改进平方根法求解结果 y:');
disp(y2_ldlt');
disp('MATLAB内置法求解结果 y:');
disp(y2_ref');

err1_chol = norm(y1_chol - y_true) / norm(y_true);
err2_chol = norm(y2_chol - y_true) / norm(y_true); 
err2_ldlt = norm(y2_ldlt - y_true) / norm(y_true);
err2_ref  = norm(y2_ref  - y_true) / norm(y_true);

res1_chol = norm(A2 * y1_chol - b2) / norm(b2);
res2_chol = norm(A2 * y2_chol - b2) / norm(b2); 
res2_ldlt = norm(A2 * y2_ldlt - b2) / norm(b2);
res2_ref  = norm(A2 * y2_ref  - b2) / norm(b2);

fprintf('真实解为全1向量 x = ones(%d,1)\n', n2);
fprintf('------------------------------------------------------------\n');
fprintf('方法类型             | 相对误差 (Error) | 相对残差 (Residual)\n');
fprintf('------------------------------------------------------------\n');
fprintf('经典平方根法         | %.4e       | %.4e\n', err1_chol, res1_chol);
fprintf('算法1.3.1平方根法    | %.4e       | %.4e\n', err2_chol, res2_chol);
fprintf('改进平方根法 (LDL'')  | %.4e       | %.4e\n', err2_ldlt, res2_ldlt);
fprintf('MATLAB 内置 \\        | %.4e       | %.4e\n', err2_ref, res2_ref);
fprintf('------------------------------------------------------------\n');

fprintf('\n====== 第(2)小题 求解时间对比 ======\n');
fprintf('经典的平方根法: %10.6f 秒\n', time_chol2);
fprintf('算法1.3.1的平方根法: %10.6f 秒\n', time_ldlt2);
fprintf('改进平方根法: %10.6f 秒\n', time_ildlt2);
fprintf('MATLAB内置法: %10.6f 秒\n', time_mat2);