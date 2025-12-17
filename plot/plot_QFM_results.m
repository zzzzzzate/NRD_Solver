% vis_QFM_stellarator.m
clear; clc; close all;

%% 1. 读取庞加莱背景数据
fprintf('Loading Poincare data...\n');
try
    poincare_data = load('poincare.dat');
    theta_back = poincare_data(:, 1);
    psi_back = poincare_data(:, 2);
    % 转换到 r/b 坐标 (r ~ sqrt(psi))
    r_back = 0.87 * sqrt(psi_back);
catch
    warning('poincare.dat not found. Run Fortran code first.');
    theta_back = []; r_back = [];
end

%% 2. 读取 QFM 曲面数据
fprintf('Loading QFM surface data...\n');
try
    qfm_data = dlmread('qfm_surfaces.dat', '', 1, 0); % 跳过第一行
    % Columns: 1:Theta, 2:Psi, 3:Nu, 4:P, 5:Q
    
    % 找到所有唯一的 (p,q) 组合
    u_pq = unique(qfm_data(:, 4:5), 'rows');
    
    % 按旋转变换 iota = p/q 的大小进行排序
    % 这对于混沌坐标系的绘制至关重要
    iotas = u_pq(:,1) ./ u_pq(:,2);
    [~, sort_idx] = sort(iotas);
    u_pq = u_pq(sort_idx, :);
    
    num_surfaces = size(u_pq, 1);
    surfaces = cell(num_surfaces, 1);
    
    fprintf('Found %d unique QFM surfaces:\n', num_surfaces);
    
    for k = 1:num_surfaces
        p = u_pq(k, 1);
        q = u_pq(k, 2);
        fprintf('  Surface %d: %d/%d (iota = %.5f)\n', k, p, q, p/q);
        
        idx = (qfm_data(:,4) == p) & (qfm_data(:,5) == q);
        dat_k = qfm_data(idx, :);
        
        % 排序 Theta 以便绘图
        [~, sort_theta_id] = sort(dat_k(:,1));
        surfaces{k} = dat_k(sort_theta_id, :);
    end
catch
    warning('qfm_surfaces.dat not found.');
    surfaces = {};
end

%% 3. 可视化配置
colors = jet(num_surfaces); % 使用 jet 颜色映射区分不同的 iota
figure('Position', [100, 100, 1400, 600], 'Color', 'w');

%% --- 子图 1: 物理空间 (Physical Space) ---
subplot(1, 2, 1);
hold on; box on; grid on;
title('Physical Space: Poincaré Plot & QFM Surfaces');
xlabel('\theta (rad)'); ylabel('r_{eff} \propto \psi^{1/2}');

% 1. 绘制庞加莱背景 (灰色点)
if ~isempty(theta_back)
    plot(theta_back, r_back,'.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 2);
end

% 2. 绘制 QFM 曲面
legend_entries = [];
legend_labels = {};

for k = 1:num_surfaces
    s = surfaces{k};
    if isempty(s), continue; end
    
    th = s(:, 1);
    psi = s(:, 2);
    r_eff = 0.87 * sqrt(psi);
    
    % 闭合曲线处理 (连接 2pi 和 0)
    th = [th; th(1) + 2*pi];
    r_eff = [r_eff; r_eff(1)];
    
    % 线型设置：根据分母 q 大小区分
    p = u_pq(k, 1); q = u_pq(k, 2);
    lw = 1.5; 
    style = '-';
    if q > 25, style = '--'; lw = 1.2; end % 高阶共振用虚线
    
    h = plot(th, r_eff, style, 'Color', colors(k,:), 'LineWidth', lw);
    
    % 仅为每种曲面添加一次图例
    legend_entries(end+1) = h;
    legend_labels{end+1} = sprintf('%d/%d (t=%.3f)', p, q, p/q);
end

xlim([0, 2*pi]);
ylim([0, 1.5]); 
legend(legend_entries, legend_labels, 'Location', 'best'); 
% 图例可能太多，建议注释掉或只显示部分

%% --- 子图 2: 混沌坐标系 (Chaotic Coordinates) ---
subplot(1, 2, 2);
hold on; box on; grid on;
title('Chaotic Coordinates: Straight Field Lines');
xlabel('\theta_{straight} (Angle)'); 
ylabel('\iota (Rotational Transform Label)');

% 在这个坐标系中，纵坐标就是 iota 值，横坐标是直线化角度
% QFM 理论保证了在这些特殊的坐标面上，磁场线是直线的。
% 图中应该显示一系列水平线。

% 1. 绘制水平线
for k = 1:num_surfaces
    s = surfaces{k};
    if isempty(s), continue; end
    
    p = u_pq(k, 1); 
    q = u_pq(k, 2);
    iota_val = p/q;
    
    % 横坐标：0 到 2pi
    th_straight = linspace(0, 2*pi, 200);
    
    % 纵坐标：恒定为 iota
    rho_straight = ones(size(th_straight)) * iota_val;
    
    style = '-';
    if q > 25, style = '--'; end
    
    plot(th_straight, rho_straight, style, 'Color', colors(k,:), 'LineWidth', 2);
    
    % 在右侧标注分数
    text(2*pi + 0.1, iota_val, sprintf('%d/%d', p, q), ...
        'Color', colors(k,:), 'FontSize', 9, 'VerticalAlignment', 'middle');
end

% 2. 示意背景混沌的"平滑化"
% (简单示意：用浅灰色填充区域)
fill([0, 2*pi, 2*pi, 0], [0.15, 0.15, 0.162, 0.162], ...
     [0.9, 0.9, 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlim([0, 2*pi]);
% 设置纵坐标范围以适应当前的 iota 范围
iota_min = min(iotas);
iota_max = max(iotas);
margin = (iota_max - iota_min) * 0.2;
if margin == 0, margin = 0.01; end
ylim([iota_min - margin, iota_max + margin]);

% 添加文本说明
text(pi, iota_min - margin/2, 'Fractal structure flattened into coordinates', ...
    'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 10);

fprintf('Visualization Complete.\n');