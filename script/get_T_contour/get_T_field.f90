module get_T_field_mod
    use precision_mod
    use boozer_mapping_params_mod
    use boozer_mapping_mod
    use sparse_solver_mod
    implicit none

contains

    ! --- Solve Heat Equation on (Psi, Theta) Grid ---
    subroutine solve_heat_equation(params, N_psi, N_theta, psi_min, psi_max, &
                                   k_par, k_perp, T_field)
        type(map_params), intent(in) :: params
        integer, intent(in) :: N_psi, N_theta
        real(dp), intent(in) :: psi_min, psi_max
        real(dp), intent(in) :: k_par, k_perp
        real(dp), intent(out) :: T_field(N_psi, N_theta)
        
        type(csr_matrix) :: A
        real(dp), allocatable :: b(:), x(:)
        real(dp) :: d_psi, d_theta, d_phi
        integer :: i, j, idx, n_dofs, istat
        real(dp) :: psi_val, theta_val, psi_f, theta_f, psi_b, theta_b
        real(dp) :: w_par, w_perp_psi, w_perp_th
        real(dp) :: err
        
        ! Arrays for constructing CSR
        integer, allocatable :: row_counts(:)
        integer, allocatable :: temp_cols(:,:)
        real(dp), allocatable :: temp_vals(:,:)
        integer, parameter :: max_entries_per_row = 40 ! 增加到40防止越界
        
        d_psi = (psi_max - psi_min) / real(N_psi - 1, dp)
        d_theta = 2.0_dp * 3.14159265358979323846_dp / real(N_theta, dp)
        d_phi = 2.0_dp * 3.14159265358979323846_dp 
        
        n_dofs = N_psi * N_theta
        
        write(*,*) 'Allocating memory for', n_dofs, 'DOFs...'
        allocate(row_counts(n_dofs), stat=istat)
        if (istat/=0) stop 'Mem alloc failed: row_counts'
        
        allocate(temp_cols(max_entries_per_row, n_dofs), stat=istat)
        if (istat/=0) stop 'Mem alloc failed: temp_cols'
        temp_cols = 0 

        allocate(temp_vals(max_entries_per_row, n_dofs), stat=istat)
        if (istat/=0) stop 'Mem alloc failed: temp_vals'
        temp_vals = 0.0_dp

        allocate(b(n_dofs), x(n_dofs), stat=istat)
        if (istat/=0) stop 'Mem alloc failed: b/x'
        
        row_counts = 0
        b = 0.0_dp
        x = 0.5_dp 
        
        w_perp_psi = k_perp / (d_psi**2)
        w_perp_th  = k_perp / (d_theta**2)
        w_par      = k_par / (d_phi**2)
        
        write(*,*) 'Building Matrix...'
        !$OMP PARALLEL DO PRIVATE(i, j, idx, psi_val, theta_val, psi_f, theta_f, psi_b, theta_b) &
        !$OMP SHARED(row_counts, temp_cols, temp_vals, b, x, params, &
        !$OMP&        N_psi, N_theta, psi_min, d_psi, d_theta, &
        !$OMP&        w_perp_psi, w_perp_th, w_par)
        do i = 1, N_psi
            do j = 1, N_theta
                idx = (i-1)*N_theta + j
                psi_val = psi_min + (i-1)*d_psi
                theta_val = (j-1)*d_theta
                
                if (i == 1) then
                    call add_entry(idx, idx, 1.0_dp, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    b(idx) = 1.0_dp
                    x(idx) = 1.0_dp
                else if (i == N_psi) then
                    call add_entry(idx, idx, 1.0_dp, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    b(idx) = 0.0_dp
                    x(idx) = 0.0_dp
                else
                    call add_entry(idx, idx, -2.0_dp*(w_perp_psi + w_perp_th) - 2.0_dp*w_par, &
                                   row_counts, temp_cols, temp_vals, max_entries_per_row)
                    
                    call add_entry(idx, idx - N_theta, w_perp_psi, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    call add_entry(idx, idx + N_theta, w_perp_psi, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    
                    if (j == 1) then
                        call add_entry(idx, idx + (N_theta-1), w_perp_th, row_counts, temp_cols, temp_vals, max_entries_per_row)
                        call add_entry(idx, idx + 1, w_perp_th, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    else if (j == N_theta) then
                        call add_entry(idx, idx - 1, w_perp_th, row_counts, temp_cols, temp_vals, max_entries_per_row)
                        call add_entry(idx, idx - (N_theta-1), w_perp_th, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    else
                        call add_entry(idx, idx - 1, w_perp_th, row_counts, temp_cols, temp_vals, max_entries_per_row)
                        call add_entry(idx, idx + 1, w_perp_th, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    end if
                    
                    call boozer_map_full(theta_val, psi_val, params, theta_f, psi_f, 1)
                    call add_interp_entries(idx, psi_f, theta_f, w_par, N_psi, N_theta, & 
                                            psi_min, d_psi, d_theta, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    
                    call boozer_map_full(theta_val, psi_val, params, theta_b, psi_b, -1)
                    call add_interp_entries(idx, psi_b, theta_b, w_par, N_psi, N_theta, & 
                                            psi_min, d_psi, d_theta, row_counts, temp_cols, temp_vals, max_entries_per_row)
                    
                    b(idx) = 0.0_dp
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        
        write(*,*) 'Assembling CSR...'
        call assemble_csr(n_dofs, row_counts, temp_cols, temp_vals, A)
        
        ! Free large temp memory
        deallocate(temp_cols, temp_vals, row_counts)
        
        write(*,*) 'Solving linear system...'
        call bicgstab(A, b, x, 1.0e-8_dp, 5000, err)
        write(*,*) 'Solver converged with relative error:', err
        
        do i = 1, N_psi
            do j = 1, N_theta
                idx = (i-1)*N_theta + j
                T_field(i, j) = x(idx)
            end do
        end do
        
        deallocate(b, x)
        
    end subroutine solve_heat_equation

    ! --- Helper: Add Matrix Entry WITH BOUNDS CHECK ---
    subroutine add_entry(row, col, val, counts, cols, vals, max_lim)
        integer, intent(in) :: row, col, max_lim
        real(dp), intent(in) :: val
        integer, intent(inout) :: counts(:)
        integer, intent(inout) :: cols(max_lim, *)
        real(dp), intent(inout) :: vals(max_lim, *)
        
        integer :: c
        
        c = counts(row) + 1
        if (c > max_lim) return ! 防止越界导致内存破坏
        
        counts(row) = c
        cols(c, row) = col
        vals(c, row) = val
    end subroutine add_entry
    
    ! --- Helper Function: Bilinear Interpolation ---
    subroutine add_interp_entries(row, psi, theta, weight, Np, Nt, p_min, dp_s, dt_s, counts, cols, vals, max_lim)
        integer, intent(in) :: row, Np, Nt, max_lim
        real(dp), intent(in) :: psi, theta, weight, p_min, dp_s, dt_s
        integer, intent(inout) :: counts(:)
        integer, intent(inout) :: cols(max_lim, *)
        real(dp), intent(inout) :: vals(max_lim, *)
        
        integer :: i0, j0, i1, j1
        real(dp) :: p_norm, t_norm, fp, ft, w00, w10, w01, w11
        integer :: idx00, idx10, idx01, idx11
        real(dp) :: th_mod
        real(dp), parameter :: TWO_PI = 6.28318530717958647692_dp
        
        ! 1. NaN 检查：如果坐标无效，直接返回，避免污染矩阵
        if (psi /= psi .or. theta /= theta) return 
        
        p_norm = (psi - p_min) / dp_s
        th_mod = modulo(theta, TWO_PI)
        t_norm = th_mod / dt_s
        
        i0 = floor(p_norm) + 1
        j0 = floor(t_norm) + 1
        
        fp = p_norm - real(i0 - 1, dp)
        ft = t_norm - real(j0 - 1, dp)
        
        ! 2. 严格的边界钳制 (Clamping)
        if (i0 < 1) then 
            i0 = 1
            fp = 0.0_dp
        elseif (i0 >= Np) then 
            i0 = Np - 1
            fp = 1.0_dp
        endif
        
        ! Theta 周期性处理 (防止 j0 越界)
        if (j0 < 1) j0 = 1
        if (j0 > Nt) j0 = Nt ! 理论上 modulo 后不会发生，但为了安全
        
        i1 = i0 + 1
        j1 = j0 + 1
        if (j1 > Nt) j1 = 1 
        
        w00 = (1.0_dp - fp) * (1.0_dp - ft)
        w10 = fp * (1.0_dp - ft)
        w01 = (1.0_dp - fp) * ft
        w11 = fp * ft
        
        ! 3. 计算索引并进行有效性检查
        idx00 = (i0-1)*Nt + j0
        idx10 = (i1-1)*Nt + j0
        idx01 = (i0-1)*Nt + j1
        idx11 = (i1-1)*Nt + j1
        
        ! 如果算出的索引非法，不要写入！
        if (idx00 < 1 .or. idx00 > Np*Nt) return
        if (idx10 < 1 .or. idx10 > Np*Nt) return
        if (idx01 < 1 .or. idx01 > Np*Nt) return
        if (idx11 < 1 .or. idx11 > Np*Nt) return
        
        call add_entry(row, idx00, weight * w00, counts, cols, vals, max_lim)
        call add_entry(row, idx10, weight * w10, counts, cols, vals, max_lim)
        call add_entry(row, idx01, weight * w01, counts, cols, vals, max_lim)
        call add_entry(row, idx11, weight * w11, counts, cols, vals, max_lim)
        
    end subroutine add_interp_entries

    subroutine assemble_csr(n, counts, cols, vals, A)
        integer, intent(in) :: n
        integer, intent(in) :: counts(:)
        integer, intent(in) :: cols(:,:)
        real(dp), intent(in) :: vals(:,:)
        type(csr_matrix), intent(out) :: A
        
        integer :: i, k
        integer(8) :: nnz_total ! Use 8-byte integer for total count check
        
        A%n = n
        nnz_total = sum(int(counts, 8))
        write(*,*) 'Total Non-zeros:', nnz_total
        A%nnz = int(nnz_total)
        
        allocate(A%row_ptr(n+1))
        allocate(A%col_ind(A%nnz))
        allocate(A%values(A%nnz))
        
        A%row_ptr(1) = 1
        do i = 1, n
            A%row_ptr(i+1) = A%row_ptr(i) + counts(i)
            do k = 1, counts(i)
                A%col_ind(A%row_ptr(i) + k - 1) = cols(k, i)
                A%values(A%row_ptr(i) + k - 1) = vals(k, i)
            end do
        end do
        
    end subroutine assemble_csr

end module get_T_field_mod