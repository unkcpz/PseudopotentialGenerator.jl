module wrapper

! wrapper functions so it can be called from Julia
! Orginal code from dftatom use modern Fortran features like allocatable arrays
! and assumed-shape arrays, which are not supported by Julia's ccall interface.

use iso_c_binding, only: c_double, c_int, c_bool
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
    integrate_adams_original => integrate_adams, &
    rk4_integrate_original => rk4_integrate, &
    get_midpoints_original => get_midpoints
use rschroed, only: &
    schroed_inward_adams_original => schroed_inward_adams, &
    schroed_outward_adams_original => schroed_outward_adams
use reigen, only: &
    solve_radial_eigenproblem_original => solve_radial_eigenproblem
use drivers, only: &
    atom_lda_original => atom_lda, &
    get_atom_orb_original => get_atom_orb
use dft, only: &
    get_Vxc_original => get_Vxc
use rpoisson, only: &
    rpoisson_outward_pc_original => rpoisson_outward_pc

implicit none

private
public mesh_exp, mesh_exp_deriv, mesh_exp_deriv2, &
    integrate_trapz_1, integrate_trapz_3, integrate_trapz_5, integrate_trapz_7, &
    integrate_simpson, integrate_adams, &
    schroed_outward_adams, schroed_inward_adams, &
    solve_radial_eigenproblem, &
    atom_lda, get_Vxc, &
    get_midpoints

contains

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

integer(c_int) function rk4_integrate(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, N) result(imax) bind(c)
    ! Wrapper for the rk4_integrate() subroutine
    integer(c_int), intent(in) :: N ! Number of grid points length(R) in Julia
    real(c_double), intent(in) :: R(N) ! Grid
    real(c_double), intent(in) :: y0(2) ! Initial condition
    ! Coefficients C1 and C2 at grid points and midpoints:
    real(c_double), intent(in) :: C1(N), C2(N), C1mid(N), C2mid(N)
    ! Maximum value (if y1 > max_val, the integration stops)
    real(c_double), intent(in) :: max_val
    ! Solution y1 and y2
    real(c_double), intent(out) :: y1(N), y2(N)

    ! The integration stops at R(imax)
    call rk4_integrate_original(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
end function

subroutine get_midpoints(R, V, N, Vmid) bind(c)
    ! Wrapper for the get_midpoints() subroutine
    integer(c_int), intent(in) :: N
    real(c_double), intent(in) :: R(N), V(N)
    real(c_double), intent(out) :: Vmid(N-1)

    Vmid = get_midpoints_original(R, V)
end subroutine

!
! Wrappers for rschroed.90
!
integer(c_int) function schroed_outward_adams(l, Z, E, R, Rp, V, P, Q, N) result(imax) bind(c)
    integer(c_int), intent(in) :: l
    integer(c_int), intent(in) :: Z
    integer(c_int), intent(in) :: N
    real(c_double), intent(in) :: E
    real(c_double), intent(in) :: R(N), Rp(N)
    real(c_double), intent(in) :: V(N)
    real(c_double), intent(out) :: P(N), Q(N)

    call schroed_outward_adams_original(l, Z, E, R, Rp, V, P, Q, imax)
end function

integer(c_int) function schroed_inward_adams(l, E, R, Rp, V, P, Q, N) result(imin) bind(c)
    integer(c_int), intent(in) :: l
    integer(c_int), intent(in) :: N
    real(c_double), intent(in) :: E
    real(c_double), intent(in) :: R(N), Rp(N)
    real(c_double), intent(in) :: V(N)
    real(c_double), intent(out) :: P(N), Q(N)

    call schroed_inward_adams_original(l, E, R, Rp, V, P, Q, imin)
end function

subroutine solve_radial_eigenproblem(n, l, Ein, eps, max_iter, &
    R, Rp, V, Z, c, relat, perturb, Emin_init, Emax_init, converged, E, P, Q, NN) bind(c)
    integer(c_int), intent(in) :: n, l, relat, Z, max_iter
    integer(c_int), intent(in) :: NN
    real(c_double), intent(in) :: R(NN), Rp(NN), V(NN), eps, Ein, c
    logical(c_bool), intent(in) :: perturb
    real(c_double), intent(in) :: Emin_init, Emax_init
    integer(c_int), intent(out) :: converged
    real(c_double), intent(out) :: P(NN), Q(NN), E
    logical :: fperturb

    fperturb = logical(perturb)
    call solve_radial_eigenproblem_original(n, l, Ein, eps, max_iter, &
        R, Rp, V, Z, c, relat, fperturb, Emin_init, Emax_init, converged, E, P, Q)
end subroutine

!
! Wrappers for dft.f90
!
subroutine get_Vxc(R, rho, relat, c, exc, Vxc, N) bind(c)
    integer(c_int), intent(in) :: N
    real(c_double), intent(in) :: R(N), rho(N)
    logical(c_bool), intent(in) :: relat
    real(c_double), intent(in) :: c
    real(c_double), intent(out) :: exc(N), Vxc(N)

    logical :: frelat

    frelat = logical(relat)
    call get_Vxc_original(R, rho, frelat, c, exc, Vxc)
end subroutine

!
! Wrappers for rpoisson.f90
!
subroutine rpoisson_outward_pc(R, Rp, rho, Vh, N) bind(c)
    integer(c_int), intent(in) :: N
    real(c_double), intent(in) :: R(N), Rp(N), rho(N)
    real(c_double), intent(out) :: Vh(N)

    Vh = rpoisson_outward_pc_original(R, Rp, rho)
end subroutine

! 
! Wrappers for drivers.f90
!
subroutine atom_lda(Z, r_min, r_max, a, N, n_orb, &
    no, lo, fo, ks_energies, E_tot, R, Rp, V_tot, density, orbitals, &
    reigen_eps, reigen_max_iter, mixing_eps, mixing_alpha, &
    mixing_max_iter, perturb) bind(c)
    integer(c_int), intent(in) :: Z
    real(c_double), intent(in) :: r_min, r_max, a
    integer(c_int), intent(in) :: N, n_orb
    integer(c_int), intent(out) :: no(n_orb), lo(n_orb)
    real(c_double), intent(out) :: fo(n_orb)
    real(c_double), intent(out) :: ks_energies(n_orb)
    real(c_double), intent(out) :: E_tot
    real(c_double), intent(out) :: R(N+1), Rp(N+1)
    real(c_double), intent(out) :: V_tot(N+1)
    real(c_double), intent(out) :: density(N+1)
    real(c_double), intent(out) :: orbitals(N+1, n_orb)
    real(c_double), intent(in) :: reigen_eps, mixing_eps, mixing_alpha
    integer(c_int), intent(in) :: mixing_max_iter, reigen_max_iter
    logical(c_bool), intent(in) :: perturb
    logical :: fperturb

    fperturb = logical(perturb)

    call atom_lda_original(Z, r_min, r_max, a, N, no, lo, fo, ks_energies, E_tot, R, Rp, &
        V_tot, density, orbitals, reigen_eps, reigen_max_iter, mixing_eps, &
        mixing_alpha, mixing_max_iter, fperturb)
end subroutine

end module
