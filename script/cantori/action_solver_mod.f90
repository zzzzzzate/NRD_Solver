! action_solver_mod.f90
module action_solver_mod
    use boozer_mapping_mod
    use boozer_mapping_params_mod

    implicit none
    
contains

    ! --- 计算一条完整轨道的 Action ---
    function calc_orbit_action(theta_start, psi_start, q_res, params) result(action)
        real(dp), intent(in) :: theta_start, psi_start
        integer, intent(in) :: q_res
        type(map_params), intent(in) :: params
        real(dp) :: action
        
        real(dp) :: th, ps, dzet, zeta
        real(dp) :: th_next, ps_next
        real(dp) :: h_val_curr, h_val_next
        real(dp) :: term_pdq, term_hdt
        integer :: step, total_steps
        type(map_params) :: p_loc
        
        action = 0.0_dp
        th = theta_start
        ps = psi_start
        zeta = 0.0_dp
        dzet = params%dzet
        total_steps = q_res * params%n_dzet_steps
        p_loc = params
        
        ! 预计算初始 Hamiltonian
        h_val_curr = get_H(th, max(0.0_dp, ps/p_loc%psi_g), zeta, p_loc)
        
        do step = 1, total_steps
            ! 计算下一步状态
            ps_next = ps; th_next = th
            call step_one_map(th_next, ps_next, zeta, p_loc)
            
            ! 计算下一步的 Hamiltonian
            h_val_next = get_H(th_next, max(0.0_dp, ps_next/p_loc%psi_g), zeta+dzet, p_loc)
            
            ! ---------------------------------------------------
            ! 【修改】 使用梯形公式 (Trapezoidal Rule) 计算积分
            ! S = \int p dq - \int H dt
            ! ---------------------------------------------------
            
            ! 1. p dq 项: 0.5 * (p_{n} + p_{n+1}) * (q_{n+1} - q_{n})
            term_pdq = 0.5_dp * (ps + ps_next) * (th_next - th)
            
            ! 2. H dt 项: 0.5 * (H_{n} + H_{n+1}) * dt
            term_hdt = 0.5_dp * (h_val_curr + h_val_next) * dzet
            
            action = action + (term_pdq - term_hdt)
            
            ! 更新状态
            th = th_next
            ps = ps_next
            h_val_curr = h_val_next ! 复用下一步的 H 作为当前的 H
            zeta = zeta + dzet
        end do
    end function calc_orbit_action

    ! --- 计算 Action Difference (Flux) ---
    function calc_action_diff(p, q, psi_seed_O, psi_seed_X, params) result(delta_W)
        integer, intent(in) :: p, q
        real(dp), intent(in) :: psi_seed_O, psi_seed_X
        type(map_params), intent(in) :: params
        real(dp) :: delta_W
        
        real(dp) :: ps_O, nu_O, ps_X, nu_X
        real(dp) :: action_O, action_X
        logical :: success_O, success_X
        type(map_params) :: local_params
        real(dp) :: pi_val = acos(-1.0_dp)
        
        local_params = params
        
        ! 1. 寻找 O 点 (Stable) - 假定对称线 theta=0
        call solve_qfm_variational(p, q, 0.0_dp, psi_seed_O, local_params, ps_O, nu_O, success_O)
        
        
        ! 2. 寻找 X 点 (Unstable) - 假定对称线 theta=pi/q
        ! 注意：X点位置取决于具体的岛结构，通常在 pi/q
        call solve_qfm_variational(p, q, pi_val/real(q,dp), psi_seed_X, local_params, ps_X, nu_X, success_X)
        
        
        if (success_O .and. success_X) then
            ! 计算 Action
            ! 必须使用各自找到的 nu 吗？
            ! Mather Action 通常定义在原始 map (nu=0)。如果轨道存在，nu应接近0。
            ! 这里我们计算 QFM 轨道本身的 Action 差作为近似。
            
            local_params%nu = nu_O
            action_O = calc_orbit_action(0.0_dp, ps_O, q, local_params)
            
            local_params%nu = nu_X
            action_X = calc_orbit_action(pi_val/real(q,dp), ps_X, q, local_params)
            
            delta_W = abs(action_X - action_O)
        else
            delta_W = 0.0_dp ! 失败返回 0
        end if
        
    end function calc_action_diff

    !  约束面积作用量极小化法
    subroutine solve_qfm_variational(p_res, q_res, theta_fixed, psi_guess_avg, &
                                     params, psi_sol_avg, nu_sol, success)
        integer, intent(in) :: p_res, q_res
        real(dp), intent(in) :: theta_fixed   ! 固定第一个点 theta(0)
        real(dp), intent(in) :: psi_guess_avg ! 平均 Psi 的初始猜测
        type(map_params), intent(in) :: params
        real(dp), intent(out) :: psi_sol_avg, nu_sol
        logical, intent(out) :: success
        
        ! 离散点定义
        integer :: k, iter
        real(dp), allocatable :: theta(:), psi_mid(:)
        real(dp), allocatable :: a_diag(:), b_diag(:), c_diag(:), rhs(:), d_theta(:)
        real(dp) :: dzet, zeta_mid, slope, th_mid, H_val
        real(dp) :: dS_dth_k, dS_dth_km1
        real(dp) :: d2S_dth2_k, d2S_dth2_km1, d2S_cross
        real(dp) :: dpsi_dth_slope
        real(dp) :: nu_curr, avg_grad, max_resid
        real(dp) :: pi_val = acos(-1.0_dp)
        
        dzet = 2.0_dp * pi_val / real(q_res, dp) ! 每一段的 Zeta 步长
        
        allocate(theta(0:q_res))
        allocate(psi_mid(1:q_res)) ! 存储第 k 段 (k-1 -> k) 的 Psi
        
        ! 用于 TDMA 的数组 (大小为 q_res - 1，因为 theta(0) 和 theta(q) 固定/关联)
        ! 变量是 theta(1) ... theta(q-1)
        allocate(a_diag(q_res-1), b_diag(q_res-1), c_diag(q_res-1), rhs(q_res-1), d_theta(q_res-1))
        
        ! 1. 初始化轨道 (直线猜测)
        theta(0) = theta_fixed
        do k = 1, q_res
            theta(k) = theta_fixed + real(k,dp)/real(q_res,dp) * 2.0_dp * pi_val * real(p_res,dp)
        end do
        
        nu_curr = 0.0_dp ! 初始猜测 nu
        success = .false.
        
        ! ---------------- Newton 迭代主循环 ----------------
        do iter = 1, 50
            
            ! A. 给定当前 theta 分布，计算每段的 Psi 和 导数信息
            do k = 1, q_res
                ! 几何信息
                slope = (theta(k) - theta(k-1)) / dzet
                th_mid = 0.5_dp * (theta(k) + theta(k-1))
                zeta_mid = (real(k,dp) - 0.5_dp) * dzet
                
                ! 物理反演: 找到 Psi 使得 dH/dPsi = slope
                psi_mid(k) = get_psi_from_slope(slope, th_mid, zeta_mid, params)
                
                ! 计算对 theta 的二阶导数信息 (Hessian 元素)
                ! d(psi)/d(slope) = 1 / (d2H/dpsi2)
                ! d(slope)/d(theta_k) = 1/dzet
                dpsi_dth_slope = 1.0_dp / (d2H_dpsi2(th_mid, max(1e-8_dp, psi_mid(k)/params%psi_g), zeta_mid, params) + 1e-12_dp)
                
                ! 离散作用量 S_k = Psi * (th_k - th_k-1) - H * dzet
                ! 我们需要 S 对 theta 的导数
                ! Gradient Component: dS_k / d_th_k = Psi + ... (由于变分原理，项会抵消，简化为 Psi)
                ! 实际上: dS_total / d_theta_k = Psi_mid(k) - Psi_mid(k+1) - nu * dzet (伪场修正)
                ! 此处我们直接组装 Hessian 和 RHS
            end do
            
            ! B. 组装线性方程组 (针对 theta_1 到 theta_q-1)
            ! 目标: (Psi_mid(k) - Psi_mid(k+1)) = nu * dzet
            ! 线性化: J * delta_theta = - Residual
            
            max_resid = 0.0_dp
            
            do k = 1, q_res - 1
                ! Residual = Gradient - nu_term
                ! Gradient_k = \partial S / \partial theta_k = Psi_{k} - Psi_{k+1}
                ! 约束方程: Psi_{k} - Psi_{k+1} - nu * dzet = 0
                rhs(k) = -( (psi_mid(k) - psi_mid(k+1)) - nu_curr * dzet )
                max_resid = max(max_resid, abs(rhs(k)))
                
                ! Hessian 元素 (三对角)
                ! d(Psi_k)/d(th_k) = dPsi/dSlope * (1/dzet)
                ! d(Psi_k+1)/d(th_k) = dPsi/dSlope * (-1/dzet)
                
                ! 上一区间的贡献 (k)
                dpsi_dth_slope = 1.0_dp / d2H_dpsi2(0.5_dp*(theta(k)+theta(k-1)), psi_mid(k)/params%psi_g, (k-0.5_dp)*dzet, params)
                a_diag(k) = -1.0_dp * (dpsi_dth_slope / dzet) ! 对应 theta(k-1)
                
                ! 当前点的自贡献 (k 和 k+1 区间都涉及 theta_k)
                ! Diag = d(Psi_k)/dth_k - d(Psi_k+1)/dth_k
                b_diag(k) = (dpsi_dth_slope / dzet)  ! 来自 Psi_k 部分
                
                ! 下一区间的贡献 (k+1)
                dpsi_dth_slope = 1.0_dp / d2H_dpsi2(0.5_dp*(theta(k+1)+theta(k)), & 
                                                    psi_mid(k+1)/params%psi_g, (k+0.5_dp)*dzet, params)
                b_diag(k) = b_diag(k) - (dpsi_dth_slope * (-1.0_dp/dzet)) ! 来自 -Psi_k+1 部分
                
                c_diag(k) = -1.0_dp * (dpsi_dth_slope * (1.0_dp/dzet)) ! 对应 theta(k+1)
            end do
            
            ! 修正边界条件: theta(0) 和 theta(q) 是固定的 (Periodic BC 已隐含在 q 的位置上)
            ! 这里是 Fixed Endpoint 变分 (theta_0 固定，theta_q = theta_0 + 2*pi*p 固定)
            ! 所以方程组只有 1..q-1 个变量。
            ! a_diag(1) 对应的项是 theta(0)，已知，移到 RHS (但在 Newton 步中 dTheta_0 = 0，所以不需要移)
            ! c_diag(q-1) 对应的项是 theta(q)，已知，dTheta_q = 0，不需要移。
            ! 仅仅需要把 a(1) 和 c(q-1) 设为 0 吗？不，TDMA 需要正确处理矩阵形状。
            ! 上面的循环里，k=1时，a_diag(1) 指向 x(0) 的系数。因为 dx(0)=0，这一项直接扔掉。
            ! k=q-1时，c_diag(q-1) 指向 x(q) 的系数。因为 dx(q)=0，这一项直接扔掉。
            
            ! 修正数组供 TDMA 使用 (丢弃越界引用)
            ! TDMA solver expect a(i) is coeff of x(i-1). For i=1, a(1) is unused (coeff of x0).
            ! c(i) is coeff of x(i+1). For i=n, c(n) is unused (coeff of x_n+1).
            
            ! C. 求解 delta_theta
            call tdma_solver(q_res-1, a_diag, b_diag, c_diag, rhs, d_theta)
            
            ! D. 更新 theta
            do k = 1, q_res - 1
                theta(k) = theta(k) + d_theta(k)
            end do
            
            ! E. 更新 Lagrange 乘子 nu
            ! Paper 方法： nu = < dS/dtheta >
            ! 实际上，要求 periodicity，即 sum(Gradient) = 0? 
            ! 简单的更新策略：直接取平均残差，或者如果我们要找“自然”的 QFM 面，
            ! 我们要求 Psi_q - Psi_1 (跨越周期) 也要匹配。
            ! 在固定端点问题中，nu 由总位移决定。
            ! nu_new = Average(Psi_k - Psi_{k+1}) / dzet
            avg_grad = 0.0_dp
            do k = 1, q_res - 1
               avg_grad = avg_grad + (psi_mid(k) - psi_mid(k+1))
            end do
            ! 还需要加上边界处的贡献？
            ! 更稳健的方法：nu 是为了使得 "Area constraint" 满足，或者简单的，
            ! 我们仅仅想让这条轨道闭合。
            ! 在 fixed endpoint setup 下，nu 其实由 theta(q)-theta(0) 这一总长度约束隐式决定了。
            ! 我们可以通过观察平均的径向力不平衡来调整 nu，或者暂且设 nu 为常数 (寻找特定 nu 的轨道)。
            ! 这里的 nu_sol 是输出。
            ! 根据文献，Update nu to satisfy nu = < partial_2 S + ... >
            nu_curr = nu_curr + sum(rhs) / (real(q_res-1,dp) * dzet) * 0.5_dp ! 简单松弛
            
            if (max_resid < 1.0e-10_dp) then
                success = .true.
                exit
            end if
        end do
        
        ! 3. 计算结果
        psi_sol_avg = sum(psi_mid) / real(q_res, dp)
        nu_sol = nu_curr
        
        deallocate(theta, psi_mid, a_diag, b_diag, c_diag, rhs, d_theta)
        
    end subroutine solve_qfm_variational

    ! ----------------------------------------------------------------------
    ! 辅助函数：三对角矩阵求解器 (Thomas Algorithm)
    ! 求解 A * x = d
    ! a: 下对角线, b: 主对角线, c: 上对角线
    ! ----------------------------------------------------------------------
    subroutine tdma_solver(n, a, b, c, d, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n), b(n), c(n), d(n)
        real(dp), intent(out) :: x(n)
        real(dp) :: cp(n), dp_arr(n), m
        integer :: i
        
        ! Forward elimination
        cp(1) = c(1) / b(1)
        dp_arr(1) = d(1) / b(1)
        
        do i = 2, n
            m = b(i) - a(i) * cp(i-1)
            cp(i) = c(i) / m
            dp_arr(i) = (d(i) - a(i) * dp_arr(i-1)) / m
        end do
        
        ! Backward substitution
        x(n) = dp_arr(n)
        do i = n-1, 1, -1
            x(i) = dp_arr(i) - cp(i) * x(i+1)
        end do
    end subroutine tdma_solver

end module action_solver_mod