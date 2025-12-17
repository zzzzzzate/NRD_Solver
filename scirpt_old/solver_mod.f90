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