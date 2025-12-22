! boozer_mapping_params_mod.f90
module boozer_mapping_params_mod
    use precision_mod
    implicit none
    
    type :: map_params
        real(dp) :: iota0 = 0.15_dp
        real(dp) :: eps0 = 0.5_dp
        real(dp) :: eps_t = 0.5_dp
        real(dp) :: eps_x = -0.31_dp
        real(dp) :: nu = 0.0_dp
        real(dp) :: psi_g = 1.0_dp
        integer  :: n_dzet_steps = 100 ! Reduced for performance in example, tune as needed
        real(dp) :: dzet
    end type map_params

contains
    ! --- Helper to get Hamiltonian perturbation ---
    ! V(psi, theta, zeta) = H - H0
    function get_V_pert(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t1, t2, t3, two, three, four
        
        two = 2.0_dp; three = 3.0_dp; four = 4.0_dp
        
        t1 = p%eps0/4.0_dp * ( (two*p%iota0 - 1.0_dp)*cos(two*theta - zeta) &
             + two*p%iota0 * cos(two*theta) ) * psi_t_norm * p%psi_g
             
        t2 = p%eps_t/6.0_dp * ( (three*p%iota0 - 1.0_dp)*cos(three*theta - zeta) &
             - three*p%iota0 * cos(three*theta) ) * (psi_t_norm)**1.5_dp * p%psi_g

        t3 = p%eps_x/8.0_dp * ( (four*p%iota0 - 1.0_dp)*cos(four*theta - zeta) &
             + four*p%iota0 * cos(four*theta) ) * (psi_t_norm)**2.0_dp * p%psi_g
             
        val = t1 + t2 + t3
    end function get_V_pert

    ! --- dH/dtheta ---
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

    ! --- dH/dpsi (normalized psi input) ---
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

    ! --- dH/dzeta (Partial derivative for B-field calc) ---
    function dH_dzeta(theta, psi_t_norm, zeta, p) result(val)
        real(dp), intent(in) :: theta, psi_t_norm, zeta
        type(map_params), intent(in) :: p
        real(dp) :: val
        real(dp) :: t1, t2, t3, two, three, four
        
        two = 2.0_dp; three = 3.0_dp; four = 4.0_dp
        
        ! Derivatives of cos(m*theta - zeta) wrt zeta is sin(m*theta - zeta)
        
        t1 = p%eps0/4.0_dp * ( (two*p%iota0 - 1.0_dp)*sin(two*theta - zeta) ) * psi_t_norm * p%psi_g
             
        t2 = p%eps_t/6.0_dp * ( (three*p%iota0 - 1.0_dp)*sin(three*theta - zeta) ) * (psi_t_norm)**1.5_dp * p%psi_g

        t3 = p%eps_x/8.0_dp * ( (four*p%iota0 - 1.0_dp)*sin(four*theta - zeta) ) * (psi_t_norm)**2.0_dp * p%psi_g
             
        val = t1 + t2 + t3
    end function dH_dzeta

end module boozer_mapping_params_mod