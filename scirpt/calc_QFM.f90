! calc_QFM.f90
! 编译命令参考: 
! gfortran -fopenmp -O3 precision_mod.f90 map_params_mod.f90 solver_mod.f90 calc_QFM.f90 -o calc_QFM
! 运行命令: ./calc_QFM

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
