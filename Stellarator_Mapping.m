clear; clc; close all;

%% ========================================================================
% 1. 设置参数
% =========================================================================
fprintf('======================================================\n');
fprintf('Boozer Map Poincaré Plot and Flux Calculator in MATLAB\n');
fprintf('======================================================\n\n');

% --- 映射和搜索参数 ---
params.iota0 = 0.15; 
params.eps0 = 0.5; 
params.eps_t = 0.5;
params.eps_x = -0.31; 
params.u_psi = 0.0;
params.n_dzet_steps = 3600;
params.dzet = 2 * pi / params.n_dzet_steps;

% --- 庞加莱截面图参数 ---
n_particles = 30;         % 用于绘图的轨道数量
n_iterations = 5000;       % 每个轨道的迭代次数
psi_t_start = 0.01;        % 轨道起始的 psi_t 最小值
r_over_b_strat = 0.87 * sqrt(psi_t_start);
psi_t_end = 1.0;          % 轨道起始的 psi_t 最大值
r_over_b_end = 0.87 * sqrt(psi_t_end);
%% ========================================================================
% 2. 计算并绘制庞加莱截面图
% =========================================================================
fprintf('Step 1: Generating Poincaré section plot...\n');
fprintf('  Number of particles: %d\n', n_particles);
fprintf('  Iterations per particle: %d\n\n', n_iterations);

% 初始化轨道
initial_psi_t = linspace(psi_t_start, psi_t_end, n_particles);
initial_theta = pi * ones(1, n_particles); % 所有轨道从 theta = pi 开始

% 存储所有点的数组
all_theta_points = zeros(n_particles, n_iterations);
all_psi_t_points = zeros(n_particles, n_iterations);

% 启动6个并行工作者

parpool(6);
% 并行计算轨道
parfor i = 1:n_particles
    fprintf('正在计算轨道 %d / %d...\n', i, n_particles);
    
    % 为每个并行工作者复制 params 结构体
    local_params = params; 
    
    % 初始化当前轨道
    theta = initial_theta(i);
    psi_t = initial_psi_t(i);
    
    temp_theta = zeros(1, n_iterations);
    temp_psi_t = zeros(1, n_iterations);
    
    % 对每个轨道进行迭代
    for j = 1:n_iterations
        % 应用一次映射 (m=1)
        x_mapped = boozer_map([theta; psi_t], 1, local_params);
        theta = x_mapped(1);
        psi_t = x_mapped(2);
        
        % 存储结果 (将 theta 限制在 [0, 2*pi] 范围内)
        temp_theta(j) = mod(theta, 2*pi);
        temp_psi_t(j) = psi_t;
    end
    all_theta_points(i, :) = temp_theta;
    all_psi_t_points(i, :) = temp_psi_t;
end
delete(gcp('nocreate'))

%% --- 绘制结果 ---
figure;
hold on;
box on;
grid on;

% 设置颜色映射
% colors = parula(n_particles); 
colors = jet(n_particles); 
colormap("jet")

% 从psi_t转换到r/b坐标
for i = 1:n_particles
    all_r_over_b_points(i, :) = 0.87 * sqrt(all_psi_t_points(i, :));
end

% 绘制
for i = 1:n_particles
    plot(all_theta_points(i, :), all_psi_t_points(i, :), '.', 'MarkerSize', 2);
end

% plot(6.1160, 0.87*sqrt(4.849916e-01),'ro');
% plot(2.5172, 0.87*sqrt(6.038993e-01),'b*');

title('Poincaré Section of the Boozer Map');
xlabel('\theta');
ylabel('r/b');
xlim([0, 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
axis tight
set(gca, 'FontSize', 12);

fprintf('Poincaré plot generated.\n');
fprintf('======================================================\n\n');

%% 函数部分
% boozer_map
function x_mapped = boozer_map(x_initial, m, params)
    % 对初始点 (theta, psi_t) 应用 m 次 Boozer 映射
    
    psi_t = x_initial(2);
    theta = x_initial(1);
    zeta = 0.0;
    
    % 总步数 = m * (2*pi / dzet)
    total_steps = round(m * params.n_dzet_steps);
    
    for i = 1:total_steps
        % 更新 zeta
        zeta = zeta + params.dzet;
        
        % 执行一步半隐式欧拉积分
        [psi_t, theta] = step_semi_implicit_euler(psi_t, theta, zeta, params.dzet, params);
    end
    
    x_mapped = [theta; psi_t];
end

% step_semi_implicit_euler
function [psi_t_final, theta_final] = step_semi_implicit_euler(psi_t, theta, zeta, dzet, params)
    % 半隐式欧拉积分方法的一步
    
    psi_t_old = psi_t;
    theta_old = theta;
    
    max_iter = 10;
    tol = 1.0e-14;
    
    % 对 psi_t_new 的初始猜测（显式欧拉）
    psi_t_norm_old = max(0.0, psi_t_old / params.psi_g);
    dH_dth = dH_dtheta(theta_old, psi_t_norm_old, zeta, params);
    psi_t_guess = psi_t_old - (dH_dth * params.psi_g - params.u_psi) * dzet;
    psi_t_new = max(0.0, psi_t_guess);

    % 迭代求解隐式部分
    for iter = 1:max_iter
        psi_t_prev_iter = psi_t_new;
        psi_t_norm_iter = max(0.0, psi_t_prev_iter / params.psi_g);
        
        dH_dth = dH_dtheta(theta_old, psi_t_norm_iter, zeta, params);
        psi_t_new = psi_t_old - (dH_dth * params.psi_g - params.u_psi) * dzet;
        
        if abs(psi_t_new - psi_t_prev_iter) < tol
            break;
        end
    end
    
    psi_t_final = max(0.0, psi_t_new);
    
    % 更新 theta
    psi_t_norm_final = max(0.0, psi_t_final / params.psi_g);
    dH_dpsi = dH_dpsi_t(theta_old, psi_t_norm_final, zeta, params);
    theta_final = theta_old + dH_dpsi * dzet;
end

% dH_dtheta
function val = dH_dtheta(theta, psi_t_norm, zeta, p)
    % 计算 dH/dtheta
    t1 = -p.eps0/4.0 * ( (2.0*p.iota0 - 1.0)*2.0 * sin(2.0*theta - zeta) ...
         + 2.0*p.iota0*2.0 * sin(2.0*theta) ) * psi_t_norm;
    t2 = -p.eps_t/6.0 * ( (3.0*p.iota0 - 1.0)*3.0 * sin(3.0*theta - zeta) ...
         - 3.0*p.iota0*3.0 * sin(3.0*theta) ) * (psi_t_norm)^1.5;
    t3 = -p.eps_x/8.0 * ( (4.0*p.iota0 - 1.0)*4.0 * sin(4.0*theta - zeta) ...
         + 4.0*p.iota0*4.0 * sin(4.0*theta) ) * (psi_t_norm)^2.0;
    val = t1 + t2 + t3;
end

% dH_dpsi_t
function val = dH_dpsi_t(theta, psi_t_norm, zeta, p)
    % 计算 dH/dpsi_t
    t1 = p.iota0 + p.eps0/4.0 * ( (2.0*p.iota0 - 1.0)*cos(2.0*theta - zeta) ...
         + 2.0*p.iota0 * cos(2.0*theta) );
    t2 = p.eps_t/6.0 * ( (3.0*p.iota0 - 1.0)*cos(3.0*theta - zeta) ...
         - 3.0*p.iota0 * cos(3.0*theta) ) * 1.5 * (psi_t_norm)^0.5;
    t3 = p.eps_x/8.0 * ( (4.0*p.iota0 - 1.0)*cos(4.0*theta - zeta) ...
         + 4.0*p.iota0 * cos(4.0*theta) ) * 2.0 * psi_t_norm;
    val = t1 + t2 + t3;
end

