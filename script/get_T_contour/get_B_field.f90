! get_B_field.f90
module get_B_field_mod
    use precision_mod
    use boozer_mapping_params_mod
    implicit none

contains

    ! --------------------------------------------------------------------------
    ! Calculate Magnetic Field B
    ! Based on B = grad(psi_t) x grad(theta) + grad(phi) x grad(psi_p)
    ! Assuming canonical coordinates where sqrt(g) = 1 (Cartesian-like for the map model)
    !
    ! Contravariant components:
    ! B^psi   = - d(psi_p)/d(theta)
    ! B^theta =   d(psi_p)/d(psi_t)
    ! B^phi   =   1  (Assuming standard canonical form)
    ! --------------------------------------------------------------------------
    subroutine get_B_field(psi_val, theta_val, phi_val, params, B_psi, B_theta, B_phi)
        real(dp), intent(in)  :: psi_val, theta_val, phi_val
        type(map_params), intent(in) :: params
        real(dp), intent(out) :: B_psi, B_theta, B_phi
        
        real(dp) :: psi_norm
        
        psi_norm = max(0.0_dp, psi_val / params%psi_g)
        
        ! B^psi = - dH/dtheta
        ! Note: The Hamiltonian H is equivalent to psi_p
        B_psi = -dH_dtheta(theta_val, psi_norm, phi_val, params)
        
        ! B^theta = dH/dpsi
        B_theta = dH_dpsi_t(theta_val, psi_norm, phi_val, params)
        
        ! B^phi = 1.0 (Standard for this Hamiltonian representation)
        B_phi = 1.0_dp
        
    end subroutine get_B_field

end module get_B_field_mod