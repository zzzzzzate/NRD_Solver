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
