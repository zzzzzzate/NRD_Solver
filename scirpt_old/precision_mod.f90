! precision_mod.f90
module precision_mod
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
end module precision_mod
