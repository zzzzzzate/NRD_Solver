! boozer_mapping_mod.f90
module boozer_mapping_mod
    use boozer_mapping_params_mod
    implicit none

contains

! --- Boozer Map (m次迭代) ---
    subroutine boozer_map_full(theta_in, psi_in, m, p_in, theta_out, psi_out)
        real(dp), intent(in) :: theta_in, psi_in
        integer, intent(in) :: m
        type(map_params), intent(in) :: p_in
        real(dp), intent(out) :: theta_out, psi_out
        
        real(dp) :: theta, psi, zeta, dzet
        real(dp) :: psi_t_old, theta_old, psi_t_new, psi_t_guess, psi_t_prev
        real(dp) :: psi_t_norm, dH_dth, dH_dpsi
        integer :: step, total_steps, iter
        real(dp), parameter :: tol = 1.0e-14_dp
        
        theta = theta_in
        psi = psi_in
        zeta = 0.0_dp
        dzet = p_in%dzet
        total_steps = m * p_in%n_dzet_steps
        
        do step = 1, total_steps
            zeta = zeta + dzet
            psi_t_old = psi
            theta_old = theta
            
            ! Guess
            psi_t_norm = max(0.0_dp, psi_t_old / p_in%psi_g)
            dH_dth = dH_dtheta(theta_old, psi_t_norm, zeta, p_in)
            psi_t_guess = psi_t_old - (dH_dth * p_in%psi_g - p_in%nu) * dzet
            psi_t_new = max(0.0_dp, psi_t_guess)
            
            ! Iteration (Simple fixed point)
            do iter = 1, 10
                psi_t_prev = psi_t_new
                psi_t_norm = max(0.0_dp, psi_t_prev / p_in%psi_g)
                dH_dth = dH_dtheta(theta_old, psi_t_norm, zeta, p_in)
                psi_t_new = psi_t_old - (dH_dth * p_in%psi_g - p_in%nu) * dzet
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

    ! --- 单步推进 (用于计算 Jacobian) ---
    subroutine step_one_map(theta, psi, zeta, params)
        real(dp), intent(inout) :: theta, psi
        real(dp), intent(in) :: zeta
        type(map_params), intent(in) :: params
        real(dp) :: psi_t_norm, dH_dth, dH_dpsi, dzet
        dzet = params%dzet
        
        ! 简化的单步逻辑 (Explicit Symplectic Euler 近似，用于差分)
        ! 注意：为了保持一致性，最好与 boozer_map_full 逻辑完全相同
        ! 但为了 Jacobian 差分的稳定性，这里用显式近似通常足够
        psi_t_norm = max(0.0_dp, psi / params%psi_g)
        dH_dth = dH_dtheta(theta, psi_t_norm, zeta + dzet, params) 
        psi = max(0.0_dp, psi - (dH_dth * params%psi_g - params%nu) * dzet)
        
        psi_t_norm = max(0.0_dp, psi / params%psi_g)
        dH_dpsi = dH_dpsi_t(theta, psi_t_norm, zeta + dzet, params)
        theta = theta + dH_dpsi * dzet
    end subroutine step_one_map

    ! --- 计算切线映射及 Trace (Greene's Residue) ---
    subroutine boozer_map_tangent(theta_in, psi_in, m, p_in, theta_out, psi_out, TraceM)
        real(dp), intent(in) :: theta_in, psi_in
        integer, intent(in) :: m
        type(map_params), intent(in) :: p_in
        real(dp), intent(out) :: theta_out, psi_out
        real(dp), intent(out) :: TraceM
        
        real(dp) :: th, ps, zeta, dzet
        real(dp) :: Mt(2,2), J_step(2,2), M_new(2,2)
        integer :: step, total_steps
        real(dp) :: th_p, ps_p, th_m, ps_m, eps
        
        th = theta_in
        ps = psi_in
        zeta = 0.0_dp
        dzet = p_in%dzet
        total_steps = m * p_in%n_dzet_steps
        
        ! 初始化单位阵
        Mt = 0.0_dp; Mt(1,1) = 1.0_dp; Mt(2,2) = 1.0_dp
        eps = 1.0e-7_dp
        
        do step = 1, total_steps
            ! 1. 计算单步 Jacobian (有限差分)
            ! d/dtheta
            th_p = th + eps; ps_p = ps; call step_one_map(th_p, ps_p, zeta, p_in)
            th_m = th - eps; ps_m = ps; call step_one_map(th_m, ps_m, zeta, p_in)
            J_step(1,1) = (th_p - th_m) / (2.0*eps)
            J_step(2,1) = (ps_p - ps_m) / (2.0*eps)
            
            ! d/dpsi
            th_p = th; ps_p = ps + eps; call step_one_map(th_p, ps_p, zeta, p_in)
            th_m = th; ps_m = ps - eps; call step_one_map(th_m, ps_m, zeta, p_in)
            J_step(1,2) = (th_p - th_m) / (2.0*eps)
            J_step(2,2) = (ps_p - ps_m) / (2.0*eps)
            
            ! 2. 推进主轨道
            call step_one_map(th, ps, zeta, p_in)
            zeta = zeta + dzet
            
            ! 3. 矩阵乘法 M_new = J * Mt
            M_new = matmul(J_step, Mt)
            Mt = M_new
        end do
        
        theta_out = th
        psi_out = ps
        TraceM = Mt(1,1) + Mt(2,2)
    end subroutine boozer_map_tangent

        ! ----------------------------------------------------------------------
    ! 辅助函数：已知 dTheta/dZeta，反解 Psi
    ! 原理：dTheta/dZeta = dH/dPsi(Psi, Theta, Zeta)
    ! 因为 H 对 Psi 近似是线性的 (H ~ iota*psi)，这个反解非常快
    ! ----------------------------------------------------------------------
    function get_psi_from_slope(slope, theta_bar, zeta_bar, params) result(psi_val)
        real(dp), intent(in) :: slope, theta_bar, zeta_bar
        type(map_params), intent(in) :: params
        real(dp) :: psi_val
        
        real(dp) :: psi_curr, f, df, psi_norm
        integer :: iter
        
        ! 初始猜测： slope ~ iota0
        ! H0 = iota0 * psi => dH0/dpsi = iota0 => psi 关系不直接
        ! dH/dpsi = iota0 + pert
        ! 粗略猜测 psi = 0.5 (归一化值)
        psi_curr = 0.5_dp * params%psi_g 
        
        do iter = 1, 10
            psi_norm = max(1.0e-8_dp, psi_curr / params%psi_g)
            
            ! f(psi) = dH/dpsi - slope = 0
            f = dH_dpsi_t(theta_bar, psi_norm, zeta_bar, params) - slope
            
            ! df/dpsi = d2H/dpsi2
            df = d2H_dpsi2(theta_bar, psi_norm, zeta_bar, params)
            
            ! 保护除零
            if (abs(df) < 1.0e-10_dp) df = 1.0e-10_dp
            
            psi_curr = psi_curr - f / df
            
            if (abs(f) < 1.0e-12_dp) exit
        end do
        psi_val = max(1.0e-8_dp, psi_curr)
    end function get_psi_from_slope

    ! ---------------------------------------------------------------
    ! 【新增】猜测函数：根据目标 Iota 反推 Psi
    ! 原理：寻找对称线 (theta=0, zeta=0) 处，局部 rotational transform 等于目标 iota 的 Psi 值。
    ! ---------------------------------------------------------------
    function guess_psi_from_iota(target_iota, params) result(psi_guess)
        real(dp), intent(in) :: target_iota
        type(map_params), intent(in) :: params
        real(dp) :: psi_guess
        
        ! 调用现有的反解函数，在对称点求解 dH/dPsi = Iota
        psi_guess = get_psi_from_slope(target_iota, 0.0_dp, 0.0_dp, params)
        
    end function guess_psi_from_iota
end module boozer_mapping_mod
