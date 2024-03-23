module wrapper

! wrapper functions so it can be called from Julia
! Orginal code from dftatom use modern Fortran features like allocatable arrays
! and assumed-shape arrays, which are not supported by Julia's ccall interface.

use types, only: dp
use mesh, only: &
    mesh_exp_original => mesh_exp, &
    mesh_exp_deriv_original => mesh_exp_deriv, &
    mesh_exp_deriv2_original => mesh_exp_deriv2

implicit none

private
public mesh_exp

contains

subroutine mesh_exp(r_min, r_max, a, N, output_mesh)
! Wrapper for the mesh_exp() function
real(dp), intent(in) :: r_min, r_max, a
integer, intent(in) :: N
real(dp), intent(out) :: output_mesh(N+1)
output_mesh = mesh_exp_original(r_min, r_max, a, N)
end subroutine

subroutine mesh_exp_deriv(r_min, r_max, a, N, output_Rp)
! Wrapper for the mesh_exp_deriv() function
real(dp), intent(in) :: r_min, r_max, a
integer, intent(in) :: N
real(dp), intent(out) :: output_Rp(N+1)
output_Rp = mesh_exp_deriv_original(r_min, r_max, a, N)
end subroutine

subroutine mesh_exp_deriv2(r_min, r_max, a, N, output_Rpp)
! Wrapper for the mesh_exp_deriv() function
real(dp), intent(in) :: r_min, r_max, a
integer, intent(in) :: N
real(dp), intent(out) :: output_Rpp(N+1)
output_Rpp = mesh_exp_deriv2_original(r_min, r_max, a, N)
end subroutine

end module