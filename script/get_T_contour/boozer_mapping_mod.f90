! boozer_mapping_mod.f90
module boozer_mapping_mod
    use boozer_mapping_params_mod
    implicit none

contains

    ! --- Simple Symplectic Map Step ---
    subroutine step_one_map(theta, psi, zeta, params, direction)
        real(dp), intent(inout) :: theta, psi
        real(dp), intent(in) :: zeta
        type(map_params), intent(in) :: params
        integer, intent(in) :: direction ! 1 for forward, -1 for backward
        
        real(dp) :: psi_t_norm, dH_dth, dH_dpsi, dzet
        dzet = params%dzet * real(direction, dp)
        
        ! Symplectic Euler
        psi_t_norm = max(0.0_dp, psi / params%psi_g)
        dH_dth = dH_dtheta(theta, psi_t_norm, zeta + dzet, params) 
        psi = max(0.0_dp, psi - (dH_dth * params%psi_g - params%nu) * dzet)
        
        psi_t_norm = max(0.0_dp, psi / params%psi_g)
        dH_dpsi = dH_dpsi_t(theta, psi_t_norm, zeta + dzet, params)
        theta = theta + dH_dpsi * dzet
    end subroutine step_one_map

    ! --- Integrate Field Line for Phi = 2*pi ---
    subroutine boozer_map_full(theta_in, psi_in, p_in, theta_out, psi_out, direction)
        real(dp), intent(in) :: theta_in, psi_in
        type(map_params), intent(in) :: p_in
        real(dp), intent(out) :: theta_out, psi_out
        integer, intent(in) :: direction ! 1 = forward (+2pi), -1 = backward (-2pi)
        
        real(dp) :: th, ps, z, dz_full
        integer :: step, nsteps
        
        th = theta_in
        ps = psi_in
        z = 0.0_dp
        nsteps = p_in%n_dzet_steps
        dz_full = 2.0_dp * 3.14159265358979323846_dp
        
        ! In the map params, dzet usually corresponds to 2pi/N or similar.
        ! Here we ensure we integrate exactly 2pi.
        ! We use the parameters stored in p_in, assuming p_in%dzet is consistent with p_in%n_dzet_steps for one period.
        
        do step = 1, nsteps
            call step_one_map(th, ps, z, p_in, direction)
            z = z + p_in%dzet * direction
        end do
        
        theta_out = th
        psi_out = ps
    end subroutine boozer_map_full

end module boozer_mapping_mod