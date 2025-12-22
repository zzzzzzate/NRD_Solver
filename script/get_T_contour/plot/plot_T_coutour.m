%% 初始化
clear; clc; close all;

filename = 'T_field.bin';

if ~isfile(filename)
    error('二进制文件 %s 未找到。请确保重新运行 Fortran 程序生成新格式数据。', filename);
end

%% 1. 读取二进制数据
fid = fopen(filename, 'r');
if fid == -1
    error('无法打开文件');
end

% 读取元数据 (必须与 Fortran 写入顺序一致)
% Fortran 的默认整数通常是 int32，实数是 double (real*8 / selected_real_kind(15,307))
N_psi   = fread(fid, 1, 'int32');
N_theta = fread(fid, 1, 'int32');
psi_min = fread(fid, 1, 'float64');
psi_max = fread(fid, 1, 'float64');

fprintf('读取到网格参数: N_psi = %d, N_theta = %d\n', N_psi, N_theta);
fprintf('Psi 范围: [%f, %f]\n', psi_min, psi_max);

% 读取温度场数据
% Fortran 写入是列优先 (Column-Major)，Matlab 读取也是列优先，所以直接读即可
% 数据总数
n_elements = N_psi * N_theta;
T_field_linear = fread(fid, n_elements, 'float64');

fclose(fid);

% 验证数据完整性
if length(T_field_linear) ~= n_elements
    error('错误：文件数据不完整。预期 %d 个点，实际读取 %d 个点。', n_elements, length(T_field_linear));
end

%% 2. 重构网格
% 这里的 T_field 已经是 [N_psi, N_theta] 的形状了（基于 Fortran 的存储顺序）
% Fortran: T_field(N_psi, N_theta) -> 内存中先变 N_psi (行索引)
% Matlab:  reshape(..., N_psi, N_theta) -> 也是先填列
T_grid = reshape(T_field_linear, N_psi, N_theta);

% 手动生成坐标网格
psi_vec = linspace(psi_min, psi_max, N_psi);
theta_vec = linspace(0, 2*pi, N_theta);

[Theta_grid, Psi_grid] = meshgrid(theta_vec, psi_vec);

%% 3. 绘制降采样 Contour 图 (数据量太大时加速显示)
figure('Name', 'Temperature Contour', 'Color', 'w');

% 如果网格过大，绘图会很卡，进行降采样显示
stride = 1; 
if N_psi > 1024 || N_theta > 1024
    stride = ceil(max(N_psi, N_theta) / 1024);
    fprintf('网格过大，执行 %dx 降采样以加速绘图...\n', stride);
end

contourf(Theta_grid(1:stride:end, 1:stride:end), ...
         Psi_grid(1:stride:end, 1:stride:end), ...
         T_grid(1:stride:end, 1:stride:end), ...
         50, 'LineColor', 'none'); 

colorbar;
colormap('jet');
xlabel('\theta (rad)', 'FontSize', 12);
ylabel('\psi', 'FontSize', 12);
title(sprintf('Temperature Contour (%dx%d)', N_psi, N_theta), 'FontSize', 14);
axis tight;

%% 4. 绘制 T 随 Psi 的图 (Psi < 1)
target_idx = psi_vec <= 1.0;

if ~any(target_idx)
    target_idx = true(size(psi_vec)); % Fallback
end

psi_sub = psi_vec(target_idx);
T_sub = T_grid(target_idx, :);

% 计算均值
T_mean = mean(T_sub, 2);

figure('Name', 'T vs Psi Profile', 'Color', 'w');
hold on; grid on;

% 绘制散点 (同样降采样，否则点太密集变成纯色块)
scatter_stride = 5; 
psi_scatter = Psi_grid(target_idx, 1:scatter_stride:end);
T_scatter   = T_grid(target_idx, 1:scatter_stride:end);

scatter(psi_scatter(:), T_scatter(:), 5, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.5);
plot(psi_sub, T_mean, 'r-', 'LineWidth', 2);

xlim([min(psi_sub), max(psi_sub)]);
xlabel('\psi');
ylabel('Temperature');
title('T vs \psi Profile');
legend({'T distribution', 'Mean T'}, 'Location', 'Best');

hold off;