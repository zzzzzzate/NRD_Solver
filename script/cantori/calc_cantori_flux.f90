program calc_cantori
    use action_solver_mod
    use boozer_mapping_mod
    use boozer_mapping_params_mod
    ! 编译示例:
    ! gfortran -o calc_cantori calc_cantori_flux.f90 action_solver_mod.f90 boozer_mapping_mod.f90 boozer_mapping_params_mod.f90 precision_mod.f90
    implicit none
    
    type(map_params) :: params
    integer :: p1, q1, p2, q2
    integer :: max_level
    real(dp) :: psi_guess_left, psi_guess_right
    real(dp) :: nu_guess_left, nu_guess_right
    real(dp) :: iota1, iota2
    
    ! 设置高精度
    params%n_dzet_steps = 3600 
    
    ! --- 初始区间选择 ---
    ! Parent 1: 3/19 
    p1 = 3; q1 = 19 
    iota1 = real(p1,dp)/real(q1,dp)
    
    ! 使用自动猜测函数来获得 psi_guess
    psi_guess_left = guess_psi_from_iota(iota1, params)
    nu_guess_left = 0.0_dp
    print *, "Auto-guess Psi for ", p1, "/", q1, " (Iota=", iota1, ") is ", psi_guess_left

    ! Parent 2: 5/32 
    p2 = 5; q2 = 32
    iota2 = real(p2,dp)/real(q2,dp)
    
    ! 【修改】使用自动猜测函数
    psi_guess_right = guess_psi_from_iota(iota2, params)
    nu_guess_right = 0.0_dp
    print *, "Auto-guess Psi for ", p2, "/", q2, " (Iota=", iota2, ") is ", psi_guess_right
    
    max_level = 6
    
    print *, "==============================================="
    print *, " Analyzing Cantori Structure via Farey Tree"
    print *, " Parents: ", p2, "/", q2, " and ", p1, "/", q1
    print *, "==============================================="
    print '(A6, A10, A12, A15, A15)', "Level", "P/Q", "Iota", "Residue", "Delta_W"
    
    ! 调用递归分析
    call analyze_farey(p2, q2, psi_guess_right, nu_guess_right, &
                       p1, q1, psi_guess_left, nu_guess_left, &
                       1, max_level)

contains

    recursive subroutine analyze_farey(p_L, q_L, psi_L, nu_L, &
                                       p_R, q_R, psi_R, nu_R, &
                                       curr_lvl, max_lvl)
        integer, intent(in) :: p_L, q_L, p_R, q_R, curr_lvl, max_lvl
        real(dp), intent(in) :: psi_L, nu_L, psi_R, nu_R
        
        integer :: p_mid, q_mid
        real(dp) :: residue, delta_w, iota_val
        real(dp) :: psi_sol, nu_sol, trace
        real(dp) :: psi_guess, nu_guess
        real(dp) :: dummy_th, dummy_psi
        logical :: success
        type(map_params) :: local_params
        
        if (curr_lvl > max_lvl) return
        
        ! 1. Farey Sum (中间节点)
        p_mid = p_L + p_R
        q_mid = q_L + q_R
        iota_val = real(p_mid,dp) / real(q_mid,dp)
        
        ! 2. Warm Start 猜测
        psi_guess = (real(q_L,dp)*psi_L + real(q_R,dp)*psi_R) / real(q_L+q_R,dp)
        nu_guess  = (real(q_L,dp)*nu_L + real(q_R,dp)*nu_R) / real(q_L+q_R,dp)
        
        local_params = params
        
        ! 3. 求解 QFM O点 (Theta = 0)
        call solve_qfm_variational(p_mid, q_mid, 0.0_dp, psi_guess, params, &
                                   psi_sol, nu_sol, success)
        
        if (success) then
            ! 4. 计算 Residue
            local_params%nu = nu_sol 
            call boozer_map_tangent(0.0_dp, psi_sol, q_mid, local_params, &
                                    dummy_th, dummy_psi, trace)
            residue = (2.0_dp - trace) / 4.0_dp
            
            ! 5. 计算 Delta W (Flux)
            delta_w = calc_action_diff(p_mid, q_mid, psi_sol, psi_sol, params)
            
        else
            residue = -999.0_dp
            delta_w = -1.0_dp
            psi_sol = psi_guess
            nu_sol = nu_guess
        end if
        
        if (success) then
            print '(I4, I5, "/", I4, F12.6, F14.6, E16.6)', &
                  curr_lvl, p_mid, q_mid, iota_val, residue, delta_w
        else
            print '(I4, I5, "/", I4, F12.6, A14, A16)', &
                  curr_lvl, p_mid, q_mid, iota_val, "  FAILED", "  FAILED"
        end if
        
        ! 6. 递归遍历
        call analyze_farey(p_L, q_L, psi_L, nu_L, &
                           p_mid, q_mid, psi_sol, nu_sol, &
                           curr_lvl+1, max_lvl)
                           
        call analyze_farey(p_mid, q_mid, psi_sol, nu_sol, &
                           p_R, q_R, psi_R, nu_R, &
                           curr_lvl+1, max_lvl)
        
    end subroutine analyze_farey

end program calc_cantori