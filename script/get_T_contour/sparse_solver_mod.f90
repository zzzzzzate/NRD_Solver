! sparse_solver_mod.f90
module sparse_solver_mod
    use precision_mod
    implicit none
    
    ! CSR Matrix Format
    type :: csr_matrix
        integer :: n                 ! Dimension
        integer :: nnz               ! Number of non-zeros
        integer, allocatable :: row_ptr(:)
        integer, allocatable :: col_ind(:)
        real(dp), allocatable :: values(:)
    end type csr_matrix

contains

    ! Matrix-Vector Multiplication: y = A * x
    ! 包含调试检查，这会降低一些速度，但能定位 SegFault 原因
    subroutine spmv(A, x, y)
        type(csr_matrix), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: i, k, n, col
        
        n = A%n
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, k, col) SCHEDULE(STATIC)
        do i = 1, n
            y(i) = 0.0_dp
            if (A%row_ptr(i+1) > A%row_ptr(i)) then
                do k = A%row_ptr(i), A%row_ptr(i+1)-1
                    col = A%col_ind(k)
                    
                    ! --- 调试检查开始 ---
                    if (col < 1 .or. col > n) then
                        print *, "CRITICAL ERROR: Invalid index in CSR!"
                        print *, "Row:", i, " Entry:", k, " Col_Index:", col
                        ! 强制停止，避免 SegFault
                        stop 
                    end if
                    ! --- 调试检查结束 ---
                    
                    y(i) = y(i) + A%values(k) * x(col)
                end do
            end if
        end do
        !$OMP END PARALLEL DO
    end subroutine spmv

    ! Bi-CGStab Solver (Standard implementation)
    subroutine bicgstab(A, b, x, tol, max_iter, err_out)
        type(csr_matrix), intent(in) :: A
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: tol
        integer, intent(in) :: max_iter
        real(dp), intent(out) :: err_out
        
        real(dp), allocatable :: r(:), r_hat(:), p(:), v(:), s(:), t(:)
        real(dp) :: alpha, omega, rho, rho_old, beta
        real(dp) :: b_norm, resid
        integer :: i, n
        
        n = A%n
        allocate(r(n), r_hat(n), p(n), v(n), s(n), t(n))
        
        call spmv(A, x, r)
        r = b - r
        r_hat = r
        
        rho_old = 1.0_dp
        alpha = 1.0_dp
        omega = 1.0_dp
        v = 0.0_dp
        p = 0.0_dp
        
        b_norm = norm2(b)
        if (b_norm < 1.0e-20_dp) b_norm = 1.0_dp
        
        do i = 1, max_iter
            rho = dot_product(r_hat, r)
            if (abs(rho) < 1.0e-20_dp) exit
            
            if (i == 1) then
                p = r
            else
                beta = (rho / rho_old) * (alpha / omega)
                p = r + beta * (p - omega * v)
            end if
            
            call spmv(A, p, v)
            
            alpha = rho / dot_product(r_hat, v)
            
            s = r - alpha * v
            if (norm2(s) < tol * b_norm) then
                x = x + alpha * p
                exit
            end if
            
            call spmv(A, s, t)
            
            omega = dot_product(t, s) / dot_product(t, t)
            
            x = x + alpha * p + omega * s
            r = s - omega * t
            
            err_out = norm2(r) / b_norm
            if (err_out < tol) exit
            
            rho_old = rho
        end do
        
        deallocate(r, r_hat, p, v, s, t)
    end subroutine bicgstab

end module sparse_solver_mod