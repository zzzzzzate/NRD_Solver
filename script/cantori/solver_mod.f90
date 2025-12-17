! solver_mod.f90
! solver_for_boozer_mapping_and_QFM
module solver_mod
    use map_params_mod
    implicit none
    
contains

    ! --- 完整的 Boozer Map (m次迭代) ---
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
            psi_t_guess = psi_t_old - (dH_dth * p_in%psi_g - p_in%u_psi) * dzet
            psi_t_new = max(0.0_dp, psi_t_guess)
            
            ! Iteration (Simple fixed point)
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

    ! --- 单步推进 (用于计算 Jacobian) ---
    subroutine step_one_map(theta, psi, zeta, params)
        real(dp), intent(inout) :: theta, psi
        real(dp), intent(in) :: zeta
        type(map_params), intent(in) :: params
        real(dp) :: psi_t_norm, dH_dth, dH_dpsi, dzet
        dzet = params%dzet
        
        ! 简化的单步逻辑 (Explicit Symplectic Euler 近似，用于差分)
        ! 注意：为了保持一致性，最好与 map_full 逻辑完全相同
        ! 但为了 Jacobian 差分的稳定性，这里用显式近似通常足够
        psi_t_norm = max(0.0_dp, psi / params%psi_g)
        dH_dth = dH_dtheta(theta, psi_t_norm, zeta + dzet, params) 
        psi = max(0.0_dp, psi - (dH_dth * params%psi_g - params%u_psi) * dzet)
        
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

    ! --- 计算一条完整轨道的 Action ---
    function calc_orbit_action(theta_start, psi_start, q_res, params) result(action)
        real(dp), intent(in) :: theta_start, psi_start
        integer, intent(in) :: q_res
        type(map_params), intent(in) :: params
        real(dp) :: action
        
        real(dp) :: th, ps, dzet, zeta, h_val
        real(dp) :: th_next, ps_next
        integer :: step, total_steps
        type(map_params) :: p_loc
        
        action = 0.0_dp
        th = theta_start
        ps = psi_start
        zeta = 0.0_dp
        dzet = params%dzet
        total_steps = q_res * params%n_dzet_steps
        p_loc = params
        
        do step = 1, total_steps
            ! 计算下一步状态
            ps_next = ps; th_next = th
            call step_one_map(th_next, ps_next, zeta, p_loc)
            
            ! 计算 Hamiltonian (在中点或起始点)
            ! L = p * dq/dt - H
            ! Discrete Action = ps_next * (th_next - th) - H * dzet
            h_val = get_H(th, max(0.0_dp, ps_next/p_loc%psi_g), zeta+dzet, p_loc)
            
            action = action + (ps_next * (th_next - th) - h_val * dzet)
            
            th = th_next
            ps = ps_next
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
        
        local_params = params
        
        ! 1. 寻找 O 点 (Stable) - 假定对称线 theta=0
        call solve_qfm_point(p, q, 0.0_dp, psi_seed_O, 0.0_dp, local_params, ps_O, nu_O, success_O)
        
        ! 2. 寻找 X 点 (Unstable) - 假定对称线 theta=pi/q
        ! 注意：X点位置取决于具体的岛结构，通常在 pi/q
        call solve_qfm_point(p, q, 3.1415926535_dp/real(q,dp), psi_seed_X, 0.0_dp, local_params, ps_X, nu_X, success_X)
        
        if (success_O .and. success_X) then
            ! 计算 Action
            ! 必须使用各自找到的 nu 吗？
            ! Mather Action 通常定义在原始 map (nu=0)。如果轨道存在，nu应接近0。
            ! 这里我们计算 QFM 轨道本身的 Action 差作为近似。
            
            local_params%u_psi = nu_O
            action_O = calc_orbit_action(0.0_dp, ps_O, q, local_params)
            
            local_params%u_psi = nu_X
            action_X = calc_orbit_action(3.1415926535_dp/real(q,dp), ps_X, q, local_params)
            
            delta_W = abs(action_X - action_O)
        else
            delta_W = 0.0_dp ! 失败返回 0
        end if
        
    end function calc_action_diff

    ! --- Newton Solver for QFM Point ---
    subroutine solve_qfm_point(p_res, q_res, theta_fixed, psi_guess, nu_guess, params, psi_sol, nu_sol, success)
        integer, intent(in) :: p_res, q_res
        real(dp), intent(in) :: theta_fixed, psi_guess, nu_guess
        type(map_params), intent(in) :: params
        real(dp), intent(out) :: psi_sol, nu_sol
        logical, intent(out) :: success
        
        real(dp) :: x(2), f(2), dx(2), J(2,2), detJ, invJ(2,2)
        real(dp) :: theta_end, psi_end, theta_p, psi_p, theta_m, psi_m
        type(map_params) :: p_local
        integer :: iter
        real(dp), parameter :: eps_fd = 1.0e-7_dp
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        x(1) = psi_guess
        x(2) = nu_guess
        success = .false.
        p_local = params
        
        do iter = 1, 50
            p_local%u_psi = x(2)
            call boozer_map_full(theta_fixed, x(1), q_res, p_local, theta_end, psi_end)
            
            f(1) = psi_end - x(1)
            f(2) = theta_end - (theta_fixed + 2.0_dp * pi * real(p_res, dp))
            
            if (sqrt(f(1)**2 + f(2)**2) < 1.0e-10_dp) then
                success = .true.
                exit
            end if
            
            ! Jacobian (Finite Difference)
            ! dF/dpsi
            call boozer_map_full(theta_fixed, x(1)+eps_fd, q_res, p_local, theta_p, psi_p)
            J(1,1) = (psi_p - psi_end) / eps_fd - 1.0_dp
            J(2,1) = (theta_p - theta_end) / eps_fd
            
            ! dF/dnu
            p_local%u_psi = x(2) + eps_fd
            call boozer_map_full(theta_fixed, x(1), q_res, p_local, theta_m, psi_m)
            J(1,2) = (psi_m - psi_end) / eps_fd
            J(2,2) = (theta_m - theta_end) / eps_fd
            
            detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
            if (abs(detJ) < 1.0e-14_dp) return
            
            invJ(1,1) =  J(2,2) / detJ; invJ(1,2) = -J(1,2) / detJ
            invJ(2,1) = -J(2,1) / detJ; invJ(2,2) =  J(1,1) / detJ
            
            dx(1) = invJ(1,1)*f(1) + invJ(1,2)*f(2)
            dx(2) = invJ(2,1)*f(1) + invJ(2,2)*f(2)
            x = x - dx
            
            if (x(1) < 1.0e-4_dp) x(1) = 1.0e-4_dp
            if (x(1) > 2.0_dp) x(1) = 2.0_dp
        end do
        psi_sol = x(1)
        nu_sol = x(2)
    end subroutine solve_qfm_point

end module solver_mod