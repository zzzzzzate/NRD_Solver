! map_params_mod.f90
! inlcude: Boozer_nonresonant_divertor_mapping 
! inlcude: derevative functions
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
    ! --- 计算 Hamiltonian 值 (用于 Action 计算) ---
    ! H = H0(psi) + V(psi, theta, zeta)
    ! 根据 dH_dpsi 和 dH_dtheta 反推
    function get_H(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t1, t2, t3, h0
        real(dp) :: two, three, four
        
        two = 2.0_dp; three = 3.0_dp; four = 4.0_dp
        
        ! 1. H0(psi) 部分: Integral of (iota0) dpsi
        ! 注意：之前的 dH_dpsi 包含扰动项，这里只取未扰动部分积分
        ! H0 = iota0 * psi
        ! 但由于 psi_t_norm 是归一化的，这里要注意量纲，假设 H 也是归一化的
        h0 = p%iota0 * psi_t_norm * p%psi_g 

        ! 2. 扰动势函数 V (对 dH_dtheta 积分并变号，或对 dH_dpsi_pert 积分)
        ! dH_dtheta = -dV/dtheta => V = - Integral(dH_dtheta) dtheta
        
        ! Term 1 (m=2):
        ! dH/dth part: -p%eps0/4 * [ (2i-1)*2*sin(2th-z) + 2i*2*sin(2th) ] * psi
        ! Integrate w.r.t theta:
        ! V1 = p%eps0/4 * [ (2i-1)*cos(2th-z) + 2i*cos(2th) ] * psi
        t1 = p%eps0/4.0_dp * ( (two*p%iota0 - 1.0_dp)*cos(two*theta - zeta) &
             + two*p%iota0 * cos(two*theta) ) * psi_t_norm * p%psi_g
             
        ! Term 2 (m=3):
        ! V2 = p%eps_t/6 * [ (3i-1)*cos(3th-z) - 3i*cos(3th) ] * psi^1.5
        t2 = p%eps_t/6.0_dp * ( (three*p%iota0 - 1.0_dp)*cos(three*theta - zeta) &
             - three*p%iota0 * cos(three*theta) ) * (psi_t_norm)**1.5_dp * p%psi_g

        ! Term 3 (m=4):
        ! V3 = p%eps_x/8 * [ (4i-1)*cos(4th-z) + 4i*cos(4th) ] * psi^2.0
        t3 = p%eps_x/8.0_dp * ( (four*p%iota0 - 1.0_dp)*cos(four*theta - zeta) &
             + four*p%iota0 * cos(four*theta) ) * (psi_t_norm)**2.0_dp * p%psi_g
             
        val = h0 + t1 + t2 + t3
    end function get_H

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

    ! --- 解析计算二阶导数 (用于 Tangent Map / Residue) ---
    ! d2H / dtheta2
    function d2H_dth2(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t1, t2, t3, two, three, four
        two=2.0_dp; three=3.0_dp; four=4.0_dp
        
        ! 对 dH_dtheta 再求一次导 (sin -> cos)
        t1 = -p%eps0/4.0_dp * ( (2.0_dp*p%iota0 - 1.0_dp)*two*two * cos(two*theta - zeta) &
             + two*p%iota0*two*two * cos(two*theta) ) * psi_t_norm
        
        t2 = -p%eps_t/6.0_dp * ( (3.0_dp*p%iota0 - 1.0_dp)*three*three * cos(three*theta - zeta) &
             - three*p%iota0*three*three * cos(three*theta) ) * (psi_t_norm)**1.5_dp
             
        t3 = -p%eps_x/8.0_dp * ( (4.0_dp*p%iota0 - 1.0_dp)*four*four * cos(four*theta - zeta) &
             + four*p%iota0*four*four * cos(four*theta) ) * (psi_t_norm)**2.0_dp
        val = t1 + t2 + t3
    end function d2H_dth2

    ! d2H / dpsi dtheta
    function d2H_dpsi_dth(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t1, t2, t3, two, three, four
        two=2.0_dp; three=3.0_dp; four=4.0_dp
        
        ! 对 dH_dtheta 求 dpsi (注意 psi_t_norm 的幂次降阶)
        ! d(psi)/dpsi = 1/psi_g
        t1 = -p%eps0/4.0_dp * ( (2.0_dp*p%iota0 - 1.0_dp)*two * sin(two*theta - zeta) &
             + two*p%iota0*two * sin(two*theta) ) * (1.0_dp/p%psi_g)
        
        t2 = -p%eps_t/6.0_dp * ( (3.0_dp*p%iota0 - 1.0_dp)*three * sin(three*theta - zeta) &
             - three*p%iota0*three * sin(three*theta) ) * 1.5_dp*(psi_t_norm)**0.5_dp * (1.0_dp/p%psi_g)

        t3 = -p%eps_x/8.0_dp * ( (4.0_dp*p%iota0 - 1.0_dp)*four * sin(four*theta - zeta) &
             + four*p%iota0*four * sin(four*theta) ) * 2.0_dp*(psi_t_norm) * (1.0_dp/p%psi_g)
             
        val = t1 + t2 + t3
    end function d2H_dpsi_dth
    
    ! d2H / dpsi2
    function d2H_dpsi2(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t2, t3, three, four
        three=3.0_dp; four=4.0_dp
        
        ! t1 is linear in psi, so 2nd derivative is 0
        
        t2 = p%eps_t/6.0_dp * ( (three*p%iota0 - 1.0_dp)*cos(three*theta - zeta) &
             - three*p%iota0 * cos(three*theta) ) * 0.75_dp*(psi_t_norm)**(-0.5_dp) * (1.0_dp/p%psi_g**2)
             
        t3 = p%eps_x/8.0_dp * ( (four*p%iota0 - 1.0_dp)*cos(four*theta - zeta) &
             + four*p%iota0 * cos(four*theta) ) * 2.0_dp * (1.0_dp/p%psi_g**2)
             
        val = t2 + t3
    end function d2H_dpsi2

end module map_params_mod
