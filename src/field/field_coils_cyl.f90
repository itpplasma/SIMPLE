module field_coils_cyl
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi, c
    use neo_biotsavart, only: coils_t, compute_vector_potential, &
                               compute_magnetic_field

    implicit none

    type(coils_t), pointer :: active_coils => null()

contains

    subroutine set_active_coils(coils)
        type(coils_t), target, intent(in) :: coils
        active_coils => coils
    end subroutine set_active_coils


    subroutine evaluate_cyl(R, phi, Z, A_cyl, B_cyl, Bmod, gradB, dA_cyl)
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: A_cyl(3), B_cyl(3), Bmod, gradB(3)
        real(dp), intent(out), optional :: dA_cyl(3,3)

        real(dp) :: xcart(3), A_cart(3), B_cart(3)
        real(dp) :: dxcart_dcyl(3,3), dcyl_dxcart(3,3)
        real(dp) :: cos_phi, sin_phi
        real(dp) :: h, A_cart_p(3), A_cart_m(3), B_cart_p(3), B_cart_m(3)
        real(dp) :: dA_cart(3,3), dB_cart(3,3)
        real(dp) :: xcart_p(3), xcart_m(3)
        integer :: i, j

        if (.not. associated(active_coils)) then
            error stop 'field_coils_cyl: active_coils not set'
        endif

        cos_phi = cos(phi)
        sin_phi = sin(phi)

        xcart(1) = R * cos_phi
        xcart(2) = R * sin_phi
        xcart(3) = Z

        A_cart = compute_vector_potential(active_coils, xcart)
        B_cart = compute_magnetic_field(active_coils, xcart)

        dxcart_dcyl(1,1) = cos_phi
        dxcart_dcyl(1,2) = -R * sin_phi
        dxcart_dcyl(1,3) = 0d0
        dxcart_dcyl(2,1) = sin_phi
        dxcart_dcyl(2,2) = R * cos_phi
        dxcart_dcyl(2,3) = 0d0
        dxcart_dcyl(3,1) = 0d0
        dxcart_dcyl(3,2) = 0d0
        dxcart_dcyl(3,3) = 1d0

        dcyl_dxcart(1,1) = cos_phi
        dcyl_dxcart(1,2) = sin_phi
        dcyl_dxcart(1,3) = 0d0
        dcyl_dxcart(2,1) = -sin_phi / R
        dcyl_dxcart(2,2) = cos_phi / R
        dcyl_dxcart(2,3) = 0d0
        dcyl_dxcart(3,1) = 0d0
        dcyl_dxcart(3,2) = 0d0
        dcyl_dxcart(3,3) = 1d0

        A_cyl(1) = A_cart(1) * cos_phi + A_cart(2) * sin_phi
        A_cyl(2) = (-A_cart(1) * sin_phi + A_cart(2) * cos_phi) * R
        A_cyl(3) = A_cart(3)

        B_cyl(1) = B_cart(1) * cos_phi + B_cart(2) * sin_phi
        B_cyl(2) = (-B_cart(1) * sin_phi + B_cart(2) * cos_phi) * R
        B_cyl(3) = B_cart(3)

        Bmod = sqrt(B_cart(1)**2 + B_cart(2)**2 + B_cart(3)**2)

        h = 1d-6 * max(R, 1d0)

        do i = 1, 3
            xcart_p = xcart
            xcart_m = xcart
            xcart_p(i) = xcart(i) + h
            xcart_m(i) = xcart(i) - h

            B_cart_p = compute_magnetic_field(active_coils, xcart_p)
            B_cart_m = compute_magnetic_field(active_coils, xcart_m)

            dB_cart(:, i) = (B_cart_p - B_cart_m) / (2d0 * h)
        enddo

        do i = 1, 3
            gradB(i) = 0d0
            do j = 1, 3
                gradB(i) = gradB(i) + (B_cart(j) / Bmod) * &
                    (dB_cart(j,1) * dxcart_dcyl(1,i) + &
                     dB_cart(j,2) * dxcart_dcyl(2,i) + &
                     dB_cart(j,3) * dxcart_dcyl(3,i))
            enddo
        enddo

        if (present(dA_cyl)) then
            do i = 1, 3
                xcart_p = xcart
                xcart_m = xcart
                xcart_p(i) = xcart(i) + h
                xcart_m(i) = xcart(i) - h

                A_cart_p = compute_vector_potential(active_coils, xcart_p)
                A_cart_m = compute_vector_potential(active_coils, xcart_m)

                dA_cart(:, i) = (A_cart_p - A_cart_m) / (2d0 * h)
            enddo

            dA_cyl = 0d0
            do i = 1, 3
                do j = 1, 3
                    dA_cyl(1,i) = dA_cyl(1,i) + &
                        (dA_cart(1,j) * cos_phi + dA_cart(2,j) * sin_phi) * &
                        dxcart_dcyl(j,i)
                    dA_cyl(3,i) = dA_cyl(3,i) + dA_cart(3,j) * dxcart_dcyl(j,i)
                enddo
            enddo

            do i = 1, 3
                dA_cyl(2,i) = 0d0
                do j = 1, 3
                    dA_cyl(2,i) = dA_cyl(2,i) + &
                        (-dA_cart(1,j) * sin_phi + dA_cart(2,j) * cos_phi) * &
                        R * dxcart_dcyl(j,i)
                enddo
            enddo

            dA_cyl(1,2) = dA_cyl(1,2) + &
                (-A_cart(1) * sin_phi + A_cart(2) * cos_phi)
            dA_cyl(2,1) = dA_cyl(2,1) + &
                (-A_cart(1) * sin_phi + A_cart(2) * cos_phi)
            dA_cyl(2,2) = dA_cyl(2,2) + &
                (-A_cart(1) * cos_phi - A_cart(2) * sin_phi) * R
        endif

    end subroutine evaluate_cyl

end module field_coils_cyl
