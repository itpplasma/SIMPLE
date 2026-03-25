module field_gvec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    use gvec_export_data, only: gvec_export_data_t, load_gvec_export_data, &
                                gvec_profile_a_theta, gvec_profile_a_phi, &
                                gvec_profile_da_theta_ds, gvec_profile_da_phi_ds, &
                                gvec_geom_r, gvec_geom_dlambda_ds, &
                                gvec_geom_dlambda_dt, gvec_geom_dlambda_dp, &
                                gvec_geom_dr_ds, gvec_geom_dr_dt, gvec_geom_dr_dp, &
                                gvec_geom_dz_ds, gvec_geom_dz_dt, gvec_geom_dz_dp
    use gvec_reference_coordinates, only: gvec_coordinate_system_t
    use spline_vmec_sub, only: compute_field_components

    implicit none

    private
    public :: gvec_field_t, create_gvec_field

    type, extends(magnetic_field_t) :: gvec_field_t
        type(gvec_export_data_t) :: data
        character(len=256) :: filename = ''
    contains
        procedure :: evaluate => gvec_evaluate
    end type gvec_field_t

contains

    subroutine create_gvec_field(filename, field)
        character(*), intent(in) :: filename
        class(gvec_field_t), allocatable, intent(out) :: field

        allocate (gvec_field_t :: field)
        field%filename = filename
        call load_gvec_export_data(filename, field%data)

        allocate (gvec_coordinate_system_t :: field%coords)
        select type (coords => field%coords)
        type is (gvec_coordinate_system_t)
            coords%data = field%data
        class default
            error stop 'create_gvec_field: allocation failure'
        end select
    end subroutine create_gvec_field

    subroutine gvec_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        class(gvec_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3)
        real(dp), intent(out) :: hcov(3)
        real(dp), intent(out) :: Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: profiles(5)
        real(dp) :: geometry(18)
        real(dp) :: a_theta
        real(dp) :: a_phi
        real(dp) :: da_theta_ds
        real(dp) :: da_phi_ds
        real(dp) :: r
        real(dp) :: dl_ds
        real(dp) :: dl_dt
        real(dp) :: dl_dp
        real(dp) :: sqg
        real(dp) :: bctr_vartheta
        real(dp) :: bctr_varphi
        real(dp) :: bcov_s
        real(dp) :: bcov_vartheta
        real(dp) :: bcov_varphi

        call self%data%evaluate_profiles(x(1), profiles)
        call self%data%evaluate_geometry(x, geometry)

        a_theta = profiles(gvec_profile_a_theta)
        a_phi = profiles(gvec_profile_a_phi)
        da_theta_ds = profiles(gvec_profile_da_theta_ds)
        da_phi_ds = profiles(gvec_profile_da_phi_ds)
        r = geometry(gvec_geom_r)
        dl_ds = geometry(gvec_geom_dlambda_ds)
        dl_dt = geometry(gvec_geom_dlambda_dt)
        dl_dp = geometry(gvec_geom_dlambda_dp)
        call compute_field_components(r, geometry(gvec_geom_dr_ds), &
                                      geometry(gvec_geom_dr_dt), &
                                      geometry(gvec_geom_dr_dp), &
                                      geometry(gvec_geom_dz_ds), &
                                      geometry(gvec_geom_dz_dt), &
                                      geometry(gvec_geom_dz_dp), da_theta_ds, &
                                      da_phi_ds, dl_ds, dl_dt, dl_dp, sqg, &
                                      bctr_vartheta, bctr_varphi, bcov_s, &
                                      bcov_vartheta, bcov_varphi)

        Acov(1) = a_theta * dl_ds
        Acov(2) = a_theta * (1.0_dp + dl_dt)
        Acov(3) = a_phi + a_theta * dl_dp

        Bmod = sqrt(bctr_vartheta * bcov_vartheta + bctr_varphi * bcov_varphi)

        hcov(1) = (bcov_s + bcov_vartheta * dl_ds) / Bmod
        hcov(2) = bcov_vartheta * (1.0_dp + dl_dt) / Bmod
        hcov(3) = (bcov_varphi + bcov_vartheta * dl_dp) / Bmod

        if (present(sqgBctr)) then
            sqgBctr = [0.0_dp, sqg * bctr_vartheta, sqg * bctr_varphi]
        end if
    end subroutine gvec_evaluate

end module field_gvec
