module wrapper

! wrapper functions so it can be called from Julia
! Orginal code from dftatom use modern Fortran features like allocatable arrays
! and assumed-shape arrays, which are not supported by Julia's ccall interface.

use iso_c_binding, only: c_double, c_int
use mesh, only: &
    mesh_exp_original => mesh_exp, &
    mesh_exp_deriv_original => mesh_exp_deriv, &
    mesh_exp_deriv2_original => mesh_exp_deriv2
use ode1d, only: &
    integrate_trapz_1_original => integrate_trapz_1, &
    integrate_trapz_3_original => integrate_trapz_3, &
    integrate_trapz_5_original => integrate_trapz_5, &
    integrate_trapz_7_original => integrate_trapz_7, &
    integrate_simpson_original => integrate_simpson, &
    integrate_adams_original => integrate_adams

implicit none

private
public mesh_exp, mesh_exp_deriv, mesh_exp_deriv2, &
    integrate_trapz_1, integrate_trapz_3, integrate_trapz_5, integrate_trapz_7, &
    integrate_simpson, integrate_adams

contains

!
! Wrappers for ode1d.90
!
real(c_double) function integrate_trapz_1(Rp, f, N) result(s) bind(c)
! Wrapper for the integrate_trapz_1() function
real(c_double), intent(in) :: Rp(N), f(N)
integer(c_int), intent(in) :: N
s = integrate_trapz_1_original(Rp, f)
end function

real(c_double) function integrate_trapz_3(Rp, f, N) result(s) bind(c)
! Wrapper for the integrate_trapz_3() function
real(c_double), intent(in) :: Rp(N), f(N)
integer(c_int), intent(in) :: N
s = integrate_trapz_3_original(Rp, f)
end function

real(c_double) function integrate_trapz_5(Rp, f, N) result(s) bind(c)
! Wrapper for the integrate_trapz_5() function
real(c_double), intent(in) :: Rp(N), f(N)
integer(c_int), intent(in) :: N
s = integrate_trapz_5_original(Rp, f)
end function

real(c_double) function integrate_trapz_7(Rp, f, N) result(s) bind(c)
! Wrapper for the integrate_trapz_7() function
real(c_double), intent(in) :: Rp(N), f(N)
integer(c_int), intent(in) :: N
s = integrate_trapz_7_original(Rp, f)
end function

real(c_double) function integrate_simpson(Rp, f, N) result(s) bind(c)
! Wrapper for the integrate_simpson() function
real(c_double), intent(in) :: Rp(N), f(N)
integer(c_int), intent(in) :: N
s = integrate_simpson_original(Rp, f)
end function

real(c_double) function integrate_adams(Rp, f, N) result(s) bind(c)
! Wrapper for the integrate_adams() function
real(c_double), intent(in) :: Rp(N), f(N)
integer(c_int), intent(in) :: N
s = integrate_adams_original(Rp, f)
end function

!
! Wrappers for mesh.90
!
subroutine mesh_exp(r_min, r_max, a, N, output_mesh) bind(c)
! Wrapper for the mesh_exp() function
real(c_double), intent(in) :: r_min, r_max, a
integer(c_int), intent(in) :: N
real(c_double), intent(out) :: output_mesh(N+1)
output_mesh = mesh_exp_original(r_min, r_max, a, N)
end subroutine

subroutine mesh_exp_deriv(r_min, r_max, a, N, output_Rp) bind(c)
! Wrapper for the mesh_exp_deriv() function
real(c_double), intent(in) :: r_min, r_max, a
integer(c_int), intent(in) :: N
real(c_double), intent(out) :: output_Rp(N+1)
output_Rp = mesh_exp_deriv_original(r_min, r_max, a, N)
end subroutine

subroutine mesh_exp_deriv2(r_min, r_max, a, N, output_Rpp) bind(c)
! Wrapper for the mesh_exp_deriv() function
real(c_double), intent(in) :: r_min, r_max, a
integer(c_int), intent(in) :: N
real(c_double), intent(out) :: output_Rpp(N+1)
output_Rpp = mesh_exp_deriv2_original(r_min, r_max, a, N)
end subroutine

end module