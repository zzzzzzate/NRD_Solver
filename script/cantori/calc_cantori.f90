program calc_cantori
    use map_params_mod
    use solver_mod
    implicit none
    
    type(map_params) :: params
    integer :: p1, q1, p2, q2
    integer :: max_level
    real(dp) :: psi_guess_left, psi_guess_right
    real(dp) :: nu_guess_left, nu_guess_right
    
    ! 设置高精度
    params%n_dzet_steps = 100 
    
    ! --- 初始区间选择 ---
    ! 必须选择已知的、稳定的父节点。
    ! 建议先运行 calc_QFM 确定两个已知好点的 Psi 值
    
    ! Parent 1: 3/19 (Iota = 0.15789)
    p1 = 3; q1 = 19 
    ! 假设我们知道 3/19 的大致 Psi (从你之前的运行结果或猜测)
    psi_guess_left = 0.55_dp ! 你需要根据实际庞加莱图调整这个值
    nu_guess_left = 0.0_dp

    ! Parent 2: 5/32 (Iota = 0.15625) 
    ! 注意：Farey Tree 通常要求 p1/q1 < p2/q2。
    ! 0.15625 < 0.15789，所以我们应该交换顺序，或者代码里自动处理
    p2 = 5; q2 = 32
    psi_guess_right = 0.45_dp ! 同样需要根据实际情况调整
    nu_guess_right = 0.0_dp
    
    max_level = 6
    
    print *, "================================================================================"
    print *, " Analyzing Cantori Structure via Farey Tree"
    print *, " Parents: ", p2, "/", q2, " and ", p1, "/", q1
    print *, "================================================================================"
    print '(A6, A10, A12, A15, A15)', "Level", "P/Q", "Iota", "Residue", "Delta_W"
    
    ! 调用递归分析，传入初始的猜测值
    ! 注意：为了简单，我们将 p2/q2 视为左节点(较小Iota)，p1/q1 视为右节点
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
        ! 利用 Farey 邻居的性质，中间节点的 Psi 和 Nu 大致是父节点的加权平均
        ! 这里简单使用算术平均或基于 q 的加权
        psi_guess = (real(q_L,dp)*psi_L + real(q_R,dp)*psi_R) / real(q_L+q_R,dp)
        nu_guess  = (real(q_L,dp)*nu_L + real(q_R,dp)*nu_R) / real(q_L+q_R,dp)
        
        local_params = params
        
        ! 3. 求解 QFM O点 (Theta = 0)
        ! 使用高质量的 psi_guess
        call solve_qfm_point(p_mid, q_mid, 0.0_dp, psi_guess, nu_guess, local_params, &
                             psi_sol, nu_sol, success)
        
        if (success) then
            ! 4. 计算 Residue
            ! 必须设置求解出的 nu 来计算切线映射
            local_params%u_psi = nu_sol 
            
            ! 调用 boozer_map_tangent 获取 Trace (需要你在 solver_mod 中实现)
            call boozer_map_tangent(0.0_dp, psi_sol, q_mid, local_params, &
                                    dummy_th, dummy_psi, trace)
            
            residue = (2.0_dp - trace) / 4.0_dp
            
            ! 5. 计算 Delta W (Flux)
            ! 传入当前找到的 O 点 psi 作为 O_seed
            ! 传入 O 点 psi 略微偏移作为 X_seed (假设 X 点 psi 接近 O 点)
            delta_w = calc_action_diff(p_mid, q_mid, psi_sol, psi_sol, params)
            
        else
            residue = -999.0_dp
            delta_w = -1.0_dp
            ! 如果当前层级失败，尝试简单的线性插值作为解继续向下递归(防止整树断裂)
            psi_sol = psi_guess
            nu_sol = nu_guess
        end if
        
        ! 输出格式
        if (success) then
            print '(I4, I5, "/", I4, F12.6, F14.6, E16.6)', &
                  curr_lvl, p_mid, q_mid, iota_val, residue, delta_w
        else
            print '(I4, I5, "/", I4, F12.6, A14, A16)', &
                  curr_lvl, p_mid, q_mid, iota_val, "  FAILED", "  FAILED"
        end if
        
        ! 6. 递归遍历
        ! 向左分支: (Left, Mid) -> Mid 成为新的 Right
        call analyze_farey(p_L, q_L, psi_L, nu_L, &
                           p_mid, q_mid, psi_sol, nu_sol, &
                           curr_lvl+1, max_lvl)
                           
        ! 向右分支: (Mid, Right) -> Mid 成为新的 Left
        call analyze_farey(p_mid, q_mid, psi_sol, nu_sol, &
                           p_R, q_R, psi_R, nu_R, &
                           curr_lvl+1, max_lvl)
        
    end subroutine analyze_farey

end program calc_cantori