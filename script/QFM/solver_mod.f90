! solver_mod.f90
! solver_for_boozer_mapping_and_QFM
module solver_mod
    use map_params_mod
    implicit none
    
contains

    ! --- 完整的 Boozer Map (m次迭代) ---
    ! 返回最终的 (theta, psi)
    subroutine boozer_map_full(theta_in, psi_in, m, p_in, theta_out, psi_out)
        real(dp), intent(in) :: theta_in, psi_in
        integer, intent(in) :: m
        type(map_params), intent(in) :: p_in
        real(dp), intent(out) :: theta_out, psi_out
        
        real(dp) :: theta, psi, zeta, dzet
        real(dp) :: psi_t_old, theta_old, psi_t_new, psi_t_guess, psi_t_prev
        real(dp) :: psi_t_norm, dH_dth, dH_dpsi
        integer :: i, step, total_steps, iter
        real(dp), parameter :: tol = 1.0e-14_dp
        
        theta = theta_in
        psi = psi_in
        zeta = 0.0_dp
        dzet = p_in%dzet
        total_steps = m * p_in%n_dzet_steps
        
        do step = 1, total_steps
            zeta = zeta + dzet
            
            ! --- Semi-implicit Euler Step Logic Inline for Speed ---
            psi_t_old = psi
            theta_old = theta
            
            ! Guess
            psi_t_norm = max(0.0_dp, psi_t_old / p_in%psi_g)
            dH_dth = dH_dtheta(theta_old, psi_t_norm, zeta, p_in)
            ! 注意: 这里的 p_in%u_psi 就是 QFM 的 nu
            psi_t_guess = psi_t_old - (dH_dth * p_in%psi_g - p_in%u_psi) * dzet
            psi_t_new = max(0.0_dp, psi_t_guess)
            
            ! Iteration
            do iter = 1, 10
                psi_t_prev = psi_t_new
                psi_t_norm = max(0.0_dp, psi_t_prev / p_in%psi_g)
                dH_dth = dH_dtheta(theta_old, psi_t_norm, zeta, p_in)
                psi_t_new = psi_t_old - (dH_dth * p_in%psi_g - p_in%u_psi) * dzet
                if (abs(psi_t_new - psi_t_prev) < tol) exit
            end do
            
            psi = max(0.0_dp, psi_t_new)
            
            psi_t_norm = max(0.0_dp, psi / p_in%psi_g)
            dH_dpsi = dH_dpsi_t(theta_old, psi_t_norm, zeta, p_in)
            theta = theta_old + dH_dpsi * dzet
        end do
        
        theta_out = theta
        psi_out = psi
        
    end subroutine boozer_map_full

    ! --- 新增：带切线映射的 Boozer Map ---
    ! 用于计算 Greene's Residue
    subroutine boozer_map_tangent(theta_in, psi_in, m, p_in, theta_out, psi_out, TraceM)
        real(dp), intent(in) :: theta_in, psi_in
        integer, intent(in) :: m
        type(map_params), intent(in) :: p_in
        real(dp), intent(out) :: theta_out, psi_out
        real(dp), intent(out) :: TraceM ! 最终切线矩阵的迹
        
        real(dp) :: th, ps, zeta, dzet
        real(dp) :: psi_t_norm, dH_dth, dH_dpsi, d2H_dth2, d2H_dpsi2, d2H_dthdpsi
        integer :: step, total_steps
        
        ! 切线矩阵 Mt = d(th_n, ps_n) / d(th_0, ps_0)
        ! 初始化为单位阵
        real(dp) :: Mt(2,2), J_step(2,2), M_new(2,2)
        
        th = theta_in
        ps = psi_in
        zeta = 0.0_dp
        dzet = p_in%dzet
        total_steps = m * p_in%n_dzet_steps
        
        Mt(1,1) = 1.0_dp; Mt(1,2) = 0.0_dp
        Mt(2,1) = 0.0_dp; Mt(2,2) = 1.0_dp
        
        do step = 1, total_steps
            zeta = zeta + dzet
            
            ! 注意：这里需要根据你的具体哈密顿量推导 Jacobian
            ! 假设是半隐式 Euler: 
            ! 1. psi_new = psi_old - (dH/dth)*dt + nu*dt
            ! 2. th_new  = th_old  + (dH/dpsi_new)*dt
            
            ! --- 1. 计算当前点的导数 ---
            psi_t_norm = max(0.0_dp, ps / p_in%psi_g) ! 简化处理，忽略 <0 的情况
            
            ! 这里需要 map_params_mod 提供二阶导数用于 Jacobian
            ! 假设 dH_dtheta 已经有了，你需要补充 d2H...
            ! 为简化代码，这里使用有限差分近似单步 Jacobian (生产环境请用解析式)
            
            call get_step_jacobian(th, ps, zeta, p_in, J_step)
            
            ! --- 2. 推进状态 ---
            ! (复制原有的推进逻辑)
            call step_one_map(th, ps, zeta, p_in) 
            
            ! --- 3. 推进切线矩阵 M_new = J_step * Mt ---
            M_new = matmul(J_step, Mt)
            Mt = M_new
        end do
        
        theta_out = th
        psi_out = ps
        TraceM = Mt(1,1) + Mt(2,2)
        
    end subroutine boozer_map_tangent

    ! --- 辅助：单步推进 ---
    subroutine step_one_map(theta, psi, zeta, params)
        real(dp), intent(inout) :: theta, psi
        real(dp), intent(in) :: zeta
        type(map_params), intent(in) :: params
        ! ... 这里填入 boozer_map_full 中循环体的逻辑 ...
        ! (为了代码复用，建议将 map_full 的循环体提取出来)
        real(dp) :: psi_t_norm, dH_dth, dH_dpsi, psi_t_new
        ! 简化的单步逻辑示例：
        psi_t_norm = max(0.0_dp, psi / params%psi_g)
        dH_dth = dH_dtheta(theta, psi_t_norm, zeta, params)
        psi = max(0.0_dp, psi - (dH_dth * params%psi_g - params%u_psi) * params%dzet)
        
        psi_t_norm = max(0.0_dp, psi / params%psi_g)
        dH_dpsi = dH_dpsi_t(theta, psi_t_norm, zeta, params)
        theta = theta + dH_dpsi * params%dzet
    end subroutine step_one_map

    ! --- 辅助：获取单步 Jacobian (有限差分版，为了通用性) ---
    subroutine get_step_jacobian(th, ps, zeta, params, J)
        real(dp), intent(in) :: th, ps, zeta
        type(map_params), intent(in) :: params
        real(dp), intent(out) :: J(2,2)
        real(dp) :: th_p, ps_p, th_m, ps_m, eps
        real(dp) :: t_tmp, p_tmp
        
        eps = 1.0e-7_dp
        
        ! d/dtheta
        t_tmp = th + eps; p_tmp = ps
        call step_one_map(t_tmp, p_tmp, zeta, params)
        th_p = t_tmp; ps_p = p_tmp
        
        t_tmp = th - eps; p_tmp = ps
        call step_one_map(t_tmp, p_tmp, zeta, params)
        th_m = t_tmp; ps_m = p_tmp
        
        J(1,1) = (th_p - th_m) / (2.0_dp*eps)
        J(2,1) = (ps_p - ps_m) / (2.0_dp*eps)
        
        ! d/dpsi
        t_tmp = th; p_tmp = ps + eps
        call step_one_map(t_tmp, p_tmp, zeta, params)
        th_p = t_tmp; ps_p = p_tmp
        
        t_tmp = th; p_tmp = ps - eps
        call step_one_map(t_tmp, p_tmp, zeta, params)
        th_m = t_tmp; ps_m = p_tmp
        
        J(1,2) = (th_p - th_m) / (2.0_dp*eps)
        J(2,2) = (ps_p - ps_m) / (2.0_dp*eps)
    end subroutine get_step_jacobian

    ! --- 计算 Mather's Difference in Action (Flux) ---
    ! 通过计算 X 点轨道和 O 点轨道在相空间围成的“代数面积”来近似
    ! Flux = Integral (Psi_X - Psi_O) dTheta
    function calc_action_diff(p, q, params) result(delta_W)
        integer, intent(in) :: p, q
        type(map_params), intent(in) :: params
        real(dp) :: delta_W
        
        real(dp) :: th_O, ps_O, nu_O, th_X, ps_X, nu_X
        logical :: success_O, success_X
        real(dp) :: th_curr_O, ps_curr_O, th_curr_X, ps_curr_X
        real(dp) :: integral_flux
        integer :: i, step
        type(map_params) :: local_params
        real(dp) :: th_old_O, th_old_X
        
        ! 1. 寻找 O 点 (Stable) - 通常在 theta=0 或 theta=pi/q
        local_params = params
        call solve_qfm_point(p, q, 0.0_dp, 0.5_dp, 0.0_dp, local_params, ps_O, nu_O, success_O)
        
        ! 2. 寻找 X 点 (Unstable) - 通常在 theta=pi/q
        call solve_qfm_point(p, q, 3.1415926_dp/real(q,dp), 0.5_dp, 0.0_dp, local_params, ps_X, nu_X, success_X)
        
        if (.not. (success_O .and. success_X)) then
            print *, "Failed to find O or X orbit for ", p, "/", q
            delta_W = -1.0_dp
            return
        end if
        
        ! 3. 沿轨道积分计算面积差
        ! Delta W = Sum [ (Psi_X_i - Psi_O_i) * dTheta ] 
        ! 注意：这只是一个近似，精确计算需要计算 Lagrangian Action S = Sum(A*dl)
        ! 但对于简单的比较，Flux 差通常对应于 Psi 的积分差。
        
        ! 更严谨的计算：计算 Action S = Sum_{i=0}^{q-1} F(x_i, x_{i+1})
        ! 对于 Boozer Map: S = Sum [ Psi_i * (Theta_{i+1} - Theta_i) - H(Psi, Theta)*dt ] 
        ! 这里简化计算：直接计算两条轨道下的“体积”差
        
        delta_W = 0.0_dp
        th_curr_O = 0.0_dp; ps_curr_O = ps_O
        th_curr_X = 3.1415926_dp/real(q,dp); ps_curr_X = ps_X
        
        ! 修正 nu 以确保是真实的周期轨道 (对于 Cantori, nu != 0)
        ! 对于计算 Delta W，我们通常关注的是 nu=0 的情况（Resonant Island）
        ! 或者 Cantori 的 Action 差。
        
        local_params%u_psi = nu_O ! 使用 O 点的 nu (通常 O 和 X 的 nu 应该非常接近)
        
        do i = 1, q
            ! 累加 Psi * dTheta (Action variable term)
            ! 这种计算方式对应于 Action S = Int p dq
            
            ! 记录旧值
            th_old_O = th_curr_O
            th_old_X = th_curr_X
            
            ! 推进一步
            call boozer_map_full(th_curr_O, ps_curr_O, 1, local_params, th_curr_O, ps_curr_O)
            call boozer_map_full(th_curr_X, ps_curr_X, 1, local_params, th_curr_X, ps_curr_X)
            
            ! 简单的梯形法则积分：(Psi_X - Psi_O) * 2*pi/q (平均)
            ! 或者更精确的 Action 差公式
            delta_W = delta_W + (ps_curr_X - ps_curr_O) ! 这是一个非常粗略的通量指标
        end do
        
        ! 归一化
        delta_W = abs(delta_W / real(q, dp))
        
    end function calc_action_diff

    ! --- Newton Solver for QFM Point ---
    ! 给定 theta_fixed, 寻找 (psi, nu) 使得 q 次迭代后:
    ! psi_final = psi
    ! theta_final = theta_fixed + 2*pi*p
    subroutine solve_qfm_point(p_res, q_res, theta_fixed, psi_guess, nu_guess, params, psi_sol, nu_sol, success)
        integer, intent(in) :: p_res, q_res
        real(dp), intent(in) :: theta_fixed, psi_guess, nu_guess
        type(map_params), intent(in) :: params
        real(dp), intent(out) :: psi_sol, nu_sol
        logical, intent(out) :: success ! 标志是否成功找到周期点
        
        real(dp) :: x(2), f(2), dx(2)
        real(dp) :: J(2,2), detJ, invJ(2,2)
        real(dp) :: theta_end, psi_end
        real(dp) :: theta_p, psi_p, theta_m, psi_m
        type(map_params) :: p_local
        integer :: iter
        real(dp), parameter :: eps_fd = 1.0e-7_dp
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        x(1) = psi_guess
        x(2) = nu_guess
        success = .false.
        p_local = params
        
        ! ----------- 牛顿迭代求解quasi-periodic条件 ---------------
        do iter = 1, 50
            p_local%u_psi = x(2) ! 设置 nu
            
            ! 1. 计算 F(x)
            call boozer_map_full(theta_fixed, x(1), q_res, p_local, theta_end, psi_end)
            
            f(1) = psi_end - x(1)
            f(2) = theta_end - (theta_fixed + 2.0_dp * pi * real(p_res, dp))
            
            if (sqrt(f(1)**2 + f(2)**2) < 1.0e-10_dp) then
                success = .true.
                exit
            end if
            
            ! 2. 计算 Jacobian (有限差分)
            ! dF/dpsi (x1)
            call boozer_map_full(theta_fixed, x(1)+eps_fd, q_res, p_local, theta_p, psi_p)
            J(1,1) = (psi_p - psi_end) / eps_fd - 1.0_dp ! d(psi_end - psi_in)/dpsi_in
            J(2,1) = (theta_p - theta_end) / eps_fd
            
            ! dF/dnu (x2)
            p_local%u_psi = x(2) + eps_fd
            call boozer_map_full(theta_fixed, x(1), q_res, p_local, theta_m, psi_m)
            J(1,2) = (psi_m - psi_end) / eps_fd
            J(2,2) = (theta_m - theta_end) / eps_fd
            
            ! 3. 求逆并更新
            detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
            if (abs(detJ) < 1.0e-14_dp) return ! Singular
            
            invJ(1,1) =  J(2,2) / detJ
            invJ(1,2) = -J(1,2) / detJ
            invJ(2,1) = -J(2,1) / detJ
            invJ(2,2) =  J(1,1) / detJ
            
            dx(1) = invJ(1,1)*f(1) + invJ(1,2)*f(2)
            dx(2) = invJ(2,1)*f(1) + invJ(2,2)*f(2)
            
            x = x - dx
            
            ! 限制 psi 防止非物理值
            if (x(1) < 1.0e-4_dp) x(1) = 1.0e-4_dp
            if (x(1) > 1.2_dp) x(1) = 1.0_dp
        end do
        ! ---------------------- 牛顿法迭代结束 ------------------------------------
        psi_sol = x(1)
        nu_sol = x(2)
        
    end subroutine solve_qfm_point

end module solver_mod