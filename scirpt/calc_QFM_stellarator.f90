! calc_QFM_stellarator.f90
! 编译命令参考: gfortran -fopenmp -O3 calc_QFM_stellarator.f90 -o calc_QFM
! 运行命令: ./calc_QFM

! ========================================================================
module precision_mod
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
end module precision_mod
! ========================================================================
module map_params_mod
    use precision_mod
    implicit none
    
    type :: map_params
        real(dp) :: iota0 = 0.15_dp
        real(dp) :: eps0 = 0.5_dp
        real(dp) :: eps_t = 0.5_dp
        real(dp) :: eps_x = -0.31_dp
        real(dp) :: u_psi = 0.0_dp  ! 这就是 QFM 理论中的 nu 参数
        real(dp) :: psi_g = 1.0_dp  ! 归一化因子
        integer  :: n_dzet_steps = 3600
        real(dp) :: dzet
    end type map_params

contains
    
    ! --- 导数计算 ---
    function dH_dtheta(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t1, t2, t3, two, three, four
        
        two = 2.0_dp; three = 3.0_dp; four = 4.0_dp
        
        t1 = -p%eps0/4.0_dp * ( (2.0_dp*p%iota0 - 1.0_dp)*two * sin(two*theta - zeta) &
             + two*p%iota0*two * sin(two*theta) ) * psi_t_norm
        
        t2 = -p%eps_t/6.0_dp * ( (3.0_dp*p%iota0 - 1.0_dp)*three * sin(three*theta - zeta) &
             - three*p%iota0*three * sin(three*theta) ) * (psi_t_norm)**1.5_dp
             
        t3 = -p%eps_x/8.0_dp * ( (4.0_dp*p%iota0 - 1.0_dp)*four * sin(four*theta - zeta) &
             + four*p%iota0*four * sin(four*theta) ) * (psi_t_norm)**2.0_dp
             
        val = t1 + t2 + t3
    end function dH_dtheta

    function dH_dpsi_t(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t1, t2, t3, two, three, four
        
        two = 2.0_dp; three = 3.0_dp; four = 4.0_dp
        
        t1 = p%iota0 + p%eps0/4.0_dp * ( (two*p%iota0 - 1.0_dp)*cos(two*theta - zeta) &
             + two*p%iota0 * cos(two*theta) )
             
        t2 = p%eps_t/6.0_dp * ( (three*p%iota0 - 1.0_dp)*cos(three*theta - zeta) &
             - three*p%iota0 * cos(three*theta) ) * 1.5_dp * (psi_t_norm)**0.5_dp
             
        t3 = p%eps_x/8.0_dp * ( (four*p%iota0 - 1.0_dp)*cos(four*theta - zeta) &
             + four*p%iota0 * cos(four*theta) ) * two * psi_t_norm
             
        val = t1 + t2 + t3
    end function dH_dpsi_t

end module map_params_mod

! ========================================================================
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
! =======================================================================================

program calc_QFM
    use map_params_mod
    use solver_mod
    use omp_lib
    implicit none
    
    type(map_params) :: params
    integer :: n_particles, n_iter, i, j, k
    real(dp) :: pi
    real(dp), allocatable :: theta0(:), psi0(:)
    real(dp) :: th_curr, psi_curr
    type(map_params) :: local_params
    
    ! QFM 变量
    integer :: n_qfm_angles = 500 ! 每个 QFM 曲面上的角度数，提高以获得更平滑的曲面
    integer, parameter :: max_qfm = 10 ! 最大支持的 QFM 曲面数
    integer, parameter :: n_surfaces = 8 ! 实际计算的 QFM 曲面数
    integer :: p_vec(max_qfm), q_vec(max_qfm)
    real(dp) :: qfm_theta, qfm_psi, qfm_nu
    real(dp) :: psi_guess_seed
    logical :: success
    
    pi = 4.0_dp * atan(1.0_dp)
    params%dzet = 2.0_dp * pi / real(params%n_dzet_steps, dp)
    
    ! =========================================================================
    ! 1. 生成庞加莱图数据 (用于matlab绘图的背景)
    ! =========================================================================
    print *, "Step 1: Generating Poincare Plot Data..."
    n_particles = 60
    n_iter = 2000
    
    open(10, file='poincare.dat', status='replace')
    
    allocate(theta0(n_particles), psi0(n_particles))
    ! 生成初始条件 (从 psi=0.01 到 1.0)
    do i = 1, n_particles
        theta0(i) = pi
        psi0(i) = 0.01_dp + (0.99_dp) * real(i-1,dp)/real(n_particles-1,dp)
    end do
    
    ! 并行计算庞加莱截面
    !$OMP PARALLEL DO PRIVATE(i, j, th_curr, psi_curr, local_params)
    do i = 1, n_particles
        th_curr = theta0(i)
        psi_curr = psi0(i)
        local_params = params
        
        ! ----------- 每个粒子跑 n_iter 圈 ---------------
        do j = 1, n_iter

            call boozer_map_full(th_curr, psi_curr, 1, local_params, th_curr, psi_curr)
            th_curr = mod(th_curr, 2.0_dp*pi)
            if (th_curr < 0.0_dp) th_curr = th_curr + 2.0_dp*pi
            
            ! 写入文件
            !$OMP CRITICAL
            write(10, *) th_curr, psi_curr
            !$OMP END CRITICAL
        end do
        ! -----------------------------------------------
        
        if (mod(i, 10) == 0) print *, "Particle ", i, " done."
    end do
    !$OMP END PARALLEL DO
    close(10)
    
    ! =========================================================================
    ! 2. QFM 曲面计算 (混沌坐标系骨架)
    ! =========================================================================
    print *, "Step 2: Calculating QFM Surfaces (Chaotic Coordinates)..."
    
    ! 目标有理数列表 (p, q) — 针对弱剪切 iota 在 [0.15,0.162]
    ! 推荐集合（覆盖从轴到边缘的窄窗）
    ! p/q: 3/20, 5/33, 2/13, 5/32, 3/19, 7/44, 4/25, 5/31
    p_vec(1) = 3;  q_vec(1) = 20  ! 0.15000
    p_vec(2) = 5;  q_vec(2) = 33  ! 0.15152
    p_vec(3) = 2;  q_vec(3) = 13  ! 0.15385
    p_vec(4) = 5;  q_vec(4) = 32  ! 0.15625
    p_vec(5) = 3;  q_vec(5) = 19  ! 0.15789 (已验证存在)
    p_vec(6) = 7;  q_vec(6) = 44  ! 0.15909
    p_vec(7) = 4;  q_vec(7) = 25  ! 0.16000
    p_vec(8) = 5;  q_vec(8) = 31  ! 0.16129
    
    open(20, file='qfm_surfaces.dat', status='replace')
    write(20, *) "Theta Psi Nu P Q"
    
    ! 对每个共振面（使用 n_surfaces）
    do k = 1, n_surfaces
        print *, "Calculating QFM for iota = ", real(p_vec(k))/real(q_vec(k))

        ! 初始 Psi 猜测：采用简单的索引比例作为稳健默认
        ! 可在失败时改为更靠谱的手动猜测或从庞加莱图读取
        psi_guess_seed = max(1.0e-4_dp, min(1.0_dp, real(k,dp) / real(n_surfaces,dp)))
        
        ! -------------------- 并行计算该面上的所有角度 ---------------------------------
        !$OMP PARALLEL DO PRIVATE(i, qfm_theta, qfm_psi, qfm_nu, success, local_params)
        do i = 1, n_qfm_angles
            qfm_theta = (real(i-1,dp)/real(n_qfm_angles,dp)) * 2.0_dp * pi
            local_params = params
            
            ! 求解 QFM 点 (Pseudo orbit starting at qfm_theta)
            call solve_qfm_point(p_vec(k), q_vec(k), qfm_theta, psi_guess_seed, 0.0_dp, &
                                 local_params, qfm_psi, qfm_nu, success)
            
            ! 如果成功找到周期轨道，则写入文件
            if (success) then
                !$OMP CRITICAL
                write(20, '(5E16.8)') qfm_theta, qfm_psi, qfm_nu, real(p_vec(k)), real(q_vec(k))
                !$OMP END CRITICAL
            end if
        end do
        !$OMP END PARALLEL DO
        ! ---------------------------------------------------------------------------
    end do
    
    close(20)
    print *, "Calculation Complete."
    
end program calc_QFM