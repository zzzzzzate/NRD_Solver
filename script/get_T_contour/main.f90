program main
    ! 并行编译示例: gfortran -fopenmp -o heat_solver main.f90 get_T_field.f90 precision_mod.f90 boozer_mapping_params_mod.f90 boozer_mapping_mod.f90 sparse_solver_mod.f90 get_B_field.f90 
    ! export OMP_NUM_THREADS=16
    ! 运行示例: ./heat_solver
    use precision_mod
    use boozer_mapping_params_mod
    use get_T_field_mod
    implicit none

    type(map_params) :: params
    integer :: N_psi, N_theta
    real(dp) :: psi_min, psi_max
    real(dp) :: k_par, k_perp
    real(dp), allocatable :: T_field(:,:)
    integer :: i, j
    
    ! --- Parameters ---
    psi_min = 0.01_dp
    psi_max = 10.0_dp
    
    ! Grid Resolution (Adjust for desired detail vs speed)
    N_psi = 512
    N_theta = 512
    
    ! Transport Coefficients
    ! k_par / k_perp ratio is critical. 
    ! Paper uses 10^10, but 10^5 is safer for initial test/convergence
    k_par = 1.0e5_dp 
    k_perp = 1.0_dp
    
    ! Mapping Params
    params%psi_g = 1.0_dp
    params%n_dzet_steps = 3600 
    params%dzet = (2.0_dp * 3.14159265358979323846_dp) / real(params%n_dzet_steps, dp)
    
    ! --- Solve ---
    allocate(T_field(N_psi, N_theta))
    
    write(*,*) 'Starting Temperature Contour Calculation...'
    write(*,*) 'Grid:', N_psi, 'x', N_theta
    write(*,*) 'Kappa Ratio:', k_par/k_perp
    
    call solve_heat_equation(params, N_psi, N_theta, psi_min, psi_max, &
                             k_par, k_perp, T_field)
                             
    ! --- Output (Revised for Binary Output) ---
    write(*,*) 'Writing binary output to T_field.bin'
    
    ! 使用 unformatted stream 写入，效率最高，且无精度丢失
    open(unit=10, file='T_field.bin', status='replace', access='stream')
    
    ! 先写入网格元数据
    write(10) N_psi
    write(10) N_theta
    write(10) psi_min
    write(10) psi_max
    
    ! 直接写入整个数组 (Fortran 列优先存储，自动处理)
    write(10) T_field
    
    close(10)
    
    deallocate(T_field)
    write(*,*) 'Done.'

end program main