module sheet_field
    !> Analytic two-volume axisymmetric torus with a current sheet at r = r0.
    !>
    !> Cylindrical (R, phi, Z), minor radius r = |(R - R0, Z)|. Units q = m = 1.
    !> The vector potential is continuous everywhere:
    !>
    !>     A_Z = -F ln R,   A_phi = psi(r)/R,   A_R = 0,
    !>
    !> with psi(r) = a r^2 + b per volume, continuous at r0 while psi' jumps.
    !> B_phi = F/R and the poloidal field (1/R) grad psi x e_phi jump in
    !> direction and slightly in magnitude across the sheet. The canonical
    !> toroidal momentum p_phi = R v_phi + psi is exactly conserved, so the
    !> full orbit reduces to two degrees of freedom z = (R, Z, p_R, p_Z) with
    !>
    !>     H = p_R^2/2 + (p_Z + F ln R)^2/2 + (p_phi - psi)^2/(2 R^2).
    !>
    !> width = 0 gives the sharp sheet: H is continuous, its gradient jumps, and
    !> the region argument selects the branch (analytically continued past the
    !> sheet). width > 0 replaces the sheet by a layer |r - r0| < width/2 in which
    !> psi blends the two branches with a C2 quintic step; outside the layer the
    !> field is the sharp one and H is smooth everywhere, so the region argument
    !> is ignored.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: field_params_t, hamiltonian, grad_hamiltonian, hess_hamiltonian
    public :: sheet_distance, grad_sheet_distance, region_of, minor_radius
    public :: bfield, psi_value

    type :: field_params_t
        real(dp) :: r_major = 10.0_dp
        real(dp) :: f_tor = 10.0_dp
        real(dp) :: r_sheet = 0.5_dp
        real(dp) :: a_in = 1.0_dp
        real(dp) :: a_out = 1.3_dp
        real(dp) :: p_phi = 0.25_dp
        real(dp) :: width = 0.0_dp
    end type field_params_t

contains

    pure real(dp) function minor_radius(par, q) result(r)
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: q(2)

        r = hypot(q(1) - par%r_major, q(2))
    end function minor_radius

    pure integer function region_of(par, q) result(region)
        !> 1 inside the sheet, 2 outside.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: q(2)

        region = merge(1, 2, minor_radius(par, q) < par%r_sheet)
    end function region_of

    pure subroutine psi_coefficients(par, region, a, b)
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(out) :: a, b

        if (region == 1) then
            a = par%a_in
            b = 0.0_dp
        else
            a = par%a_out
            b = (par%a_in - par%a_out)*par%r_sheet**2
        end if
    end subroutine psi_coefficients

    pure logical function in_layer(par, r)
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: r

        in_layer = par%width > 0.0_dp .and. &
                   abs(r - par%r_sheet) < 0.5_dp*par%width
    end function in_layer

    pure subroutine psi_radial(par, r, psi, dpsi, d2psi)
        !> psi and its r-derivatives inside the layer: the inner branch plus
        !> (outer - inner) times the quintic step S(x), x = (r - r0)/w + 1/2.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: r
        real(dp), intent(out) :: psi, dpsi, d2psi

        real(dp) :: w, x, s, s1, s2, da, d, d1, d2

        w = par%width
        x = (r - par%r_sheet)/w + 0.5_dp
        s = x**3*(10.0_dp - 15.0_dp*x + 6.0_dp*x*x)
        s1 = 30.0_dp*x*x*(x - 1.0_dp)**2/w
        s2 = 60.0_dp*x*(x - 1.0_dp)*(2.0_dp*x - 1.0_dp)/(w*w)
        da = par%a_out - par%a_in
        d = da*(r*r - par%r_sheet**2)
        d1 = 2.0_dp*da*r
        d2 = 2.0_dp*da
        psi = par%a_in*r*r + d*s
        dpsi = 2.0_dp*par%a_in*r + d1*s + d*s1
        d2psi = 2.0_dp*par%a_in + d2*s + 2.0_dp*d1*s1 + d*s2
    end subroutine psi_radial

    pure real(dp) function psi_value(par, region, q) result(psi)
        !> psi at q = (R, Z); region selects the sharp branch (width = 0).
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: q(2)

        real(dp) :: r, a, b, d1, d2

        r = minor_radius(par, q)
        if (in_layer(par, r)) then
            call psi_radial(par, r, psi, d1, d2)
        else
            call psi_coefficients(par, branch(par, region, q), a, b)
            psi = a*r*r + b
        end if
    end function psi_value

    pure integer function branch(par, region, q)
        !> Sharp field: the requested branch. Smooth field: the geometric one.
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: q(2)

        branch = region
        if (par%width > 0.0_dp) branch = region_of(par, q)
    end function branch

    pure subroutine psi_derivatives(par, region, q, psi, g, hs)
        !> psi, its gradient and Hessian in (R, Z).
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: q(2)
        real(dp), intent(out) :: psi, g(2), hs(2, 2)

        real(dp) :: x, z, r, a, b, d1, d2

        x = q(1) - par%r_major
        z = q(2)
        r = hypot(x, z)
        if (in_layer(par, r)) then
            call psi_radial(par, r, psi, d1, d2)
            g = d1*[x, z]/r
            hs(1, 1) = d2*x*x/r**2 + d1*(1.0_dp/r - x*x/r**3)
            hs(2, 2) = d2*z*z/r**2 + d1*(1.0_dp/r - z*z/r**3)
            hs(1, 2) = d2*x*z/r**2 - d1*x*z/r**3
        else
            call psi_coefficients(par, branch(par, region, q), a, b)
            psi = a*(x*x + z*z) + b
            g = 2.0_dp*a*[x, z]
            hs(1, 1) = 2.0_dp*a
            hs(2, 2) = 2.0_dp*a
            hs(1, 2) = 0.0_dp
        end if
        hs(2, 1) = hs(1, 2)
    end subroutine psi_derivatives

    pure real(dp) function hamiltonian(par, region, z) result(h)
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z(4)

        real(dp) :: u, w

        u = par%p_phi - psi_value(par, region, z(1:2))
        w = z(4) + par%f_tor*log(z(1))
        h = 0.5_dp*z(3)**2 + 0.5_dp*w*w + 0.5_dp*u*u/z(1)**2
    end function hamiltonian

    pure subroutine grad_hamiltonian(par, region, z, g)
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z(4)
        real(dp), intent(out) :: g(4)

        real(dp) :: rr, u, w, psi, gp(2), hp(2, 2)

        call psi_derivatives(par, region, z(1:2), psi, gp, hp)
        rr = z(1)
        u = par%p_phi - psi
        w = z(4) + par%f_tor*log(rr)
        g(1) = w*par%f_tor/rr - u*gp(1)/rr**2 - u*u/rr**3
        g(2) = -u*gp(2)/rr**2
        g(3) = z(3)
        g(4) = w
    end subroutine grad_hamiltonian

    pure subroutine hess_hamiltonian(par, region, z, hs)
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z(4)
        real(dp), intent(out) :: hs(4, 4)

        real(dp) :: rr, u, w, f, psi, gp(2), hp(2, 2), u_r, u_z

        call psi_derivatives(par, region, z(1:2), psi, gp, hp)
        rr = z(1)
        f = par%f_tor
        u = par%p_phi - psi
        u_r = -gp(1)
        u_z = -gp(2)
        w = z(4) + f*log(rr)
        hs = 0.0_dp
        hs(1, 1) = (f/rr)**2 - w*f/rr**2 + (u_r*u_r - u*hp(1, 1))/rr**2 &
                   - 4.0_dp*u*u_r/rr**3 + 3.0_dp*u*u/rr**4
        hs(1, 2) = (u_z*u_r - u*hp(1, 2))/rr**2 - 2.0_dp*u*u_z/rr**3
        hs(2, 1) = hs(1, 2)
        hs(2, 2) = (u_z*u_z - u*hp(2, 2))/rr**2
        hs(1, 4) = f/rr
        hs(4, 1) = hs(1, 4)
        hs(3, 3) = 1.0_dp
        hs(4, 4) = 1.0_dp
    end subroutine hess_hamiltonian

    pure real(dp) function sheet_distance(par, q) result(s)
        !> s = r - r0, the level-set function of the sheet.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: q(2)

        s = minor_radius(par, q) - par%r_sheet
    end function sheet_distance

    pure subroutine grad_sheet_distance(par, q, g)
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: q(2)
        real(dp), intent(out) :: g(2)

        real(dp) :: r

        r = minor_radius(par, q)
        g = [q(1) - par%r_major, q(2)]/r
    end subroutine grad_sheet_distance

    pure subroutine bfield(par, q, b)
        !> Physical (B_R, B_phi, B_Z) at q = (R, Z), from the local region.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: q(2)
        real(dp), intent(out) :: b(3)

        real(dp) :: psi, gp(2), hp(2, 2)

        call psi_derivatives(par, region_of(par, q), q, psi, gp, hp)
        b(1) = -gp(2)/q(1)
        b(2) = par%f_tor/q(1)
        b(3) = gp(1)/q(1)
    end subroutine bfield

end module sheet_field
