program test_canonical_full_regression
  ! Full regression test for get_canonical_coordinates
  ! Tests against known output values captured from original code
  
  use get_can_sub
  use canonical_coordinates_mod
  use new_vmec_stuff_mod
  use vector_potentail_mod
  
  ! Import nh_stencil from wherever it's defined
  integer :: nh_stencil = 3
  
  implicit none
  
  call test_canonical_full_output()
  
  print *, "All canonical full regression tests passed!"
  
contains

  subroutine test_canonical_full_output()
    ! Test that get_canonical_coordinates produces exact expected output
    double precision, parameter :: tol = 1.0d-10
    double precision :: test_value
    logical :: fullset
    double precision :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
    double precision :: sqg_test, dsqg_dr, dsqg_dt, dsqg_dp
    double precision :: B_vartheta_test, dB_vartheta_dr, dB_vartheta_dt, dB_vartheta_dp
    double precision :: B_varphi_test, dB_varphi_dr, dB_varphi_dt, dB_varphi_dp
    double precision :: G_test
    double precision :: d2sqg_rr, d2sqg_rt, d2sqg_rp, d2sqg_tt, d2sqg_tp, d2sqg_pp
    double precision :: d2bth_rr, d2bth_rt, d2bth_rp, d2bth_tt, d2bth_tp, d2bth_pp
    double precision :: d2bph_rr, d2bph_rt, d2bph_rp, d2bph_tt, d2bph_tp, d2bph_pp
    
    ! Set up configuration matching the debug run (without relying on VMEC init)
    ns = 50
    n_theta = 36
    n_phi = 81
    h_theta = 0.17951958020513104d0
    h_phi = 3.9269908169872414d-2
    hs = 2.0408163265306117d-2
    ns_s = 5
    ns_tp = 5
    
    ! Initialize torflux (needed for vector potential)
    torflux = 1.0d0  
    
    ! Allocate and initialize vector potential spline array
    if (.not. allocated(sA_phi)) then
      allocate(sA_phi(10,ns))  
      sA_phi = 0.0d0  
    end if
    
    ! Call get_canonical_coordinates
    call get_canonical_coordinates()
    
    ! Test 1: Check module variables were set correctly
    if (ns_c /= 50) then
      print *, "ERROR: ns_c incorrect. Expected: 50, Got:", ns_c
      error stop 1
    end if
    
    if (n_theta_c /= 36) then
      print *, "ERROR: n_theta_c incorrect. Expected: 36, Got:", n_theta_c
      error stop 1
    end if
    
    if (n_phi_c /= 81) then
      print *, "ERROR: n_phi_c incorrect. Expected: 81, Got:", n_phi_c
      error stop 1
    end if
    
    ! Test 2: Check spline allocation
    if (.not. allocated(s_sqg_Bt_Bp)) then
      print *, "ERROR: s_sqg_Bt_Bp not allocated"
      error stop 1
    end if
    
    if (.not. allocated(s_G_c)) then
      print *, "ERROR: s_G_c not allocated"
      error stop 1
    end if
    
    ! Test 3: Check specific spline values (from debug output)
    ! Note: These are the values after spline_can_coord completes
    
    ! Check s_sqg_Bt_Bp(1,1,1,1,1,1,1) = -17652693.862207368
    test_value = s_sqg_Bt_Bp(1,1,1,1,1,1,1)
    if (abs(test_value - (-17652693.862207368d0)) > abs(test_value) * tol) then
      print *, "ERROR: s_sqg_Bt_Bp(1,1,1,1,1,1,1) incorrect"
      print *, "Expected:", -17652693.862207368d0, "Got:", test_value
      error stop 1
    end if
    
    ! Check s_sqg_Bt_Bp(2,1,1,1,2,1,1) = 2.1729250456264708
    test_value = s_sqg_Bt_Bp(2,1,1,1,2,1,1)
    if (abs(test_value - 2.1729250456264708d0) > abs(test_value) * tol) then
      print *, "ERROR: s_sqg_Bt_Bp(2,1,1,1,2,1,1) incorrect"
      print *, "Expected:", 2.1729250456264708d0, "Got:", test_value
      error stop 1
    end if
    
    ! Check s_sqg_Bt_Bp(3,1,1,1,2,1,1) = 76691088.956212446
    test_value = s_sqg_Bt_Bp(3,1,1,1,2,1,1)
    if (abs(test_value - 76691088.956212446d0) > abs(test_value) * tol) then
      print *, "ERROR: s_sqg_Bt_Bp(3,1,1,1,2,1,1) incorrect"
      print *, "Expected:", 76691088.956212446d0, "Got:", test_value
      error stop 1
    end if
    
    ! Check s_G_c(1,1,1,1,1,1) = 1.0e-8
    test_value = s_G_c(1,1,1,1,1,1)
    if (abs(test_value - 1.0d-8) > tol) then
      print *, "ERROR: s_G_c(1,1,1,1,1,1) incorrect"
      print *, "Expected:", 1.0d-8, "Got:", test_value
      error stop 1
    end if
    
    ! Test 4: Test interpolation at a specific point
    fullset = .true.
    call splint_can_coord(fullset, 0, 0.5d0, 0.1d0, 0.1d0, &
                          A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3, &
                          sqg_test, dsqg_dr, dsqg_dt, dsqg_dp, &
                          B_vartheta_test, dB_vartheta_dr, dB_vartheta_dt, dB_vartheta_dp, &
                          B_varphi_test, dB_varphi_dr, dB_varphi_dt, dB_varphi_dp, &
                          d2sqg_rr, d2sqg_rt, d2sqg_rp, d2sqg_tt, d2sqg_tp, d2sqg_pp, &
                          d2bth_rr, d2bth_rt, d2bth_rp, d2bth_tt, d2bth_tp, d2bth_pp, &
                          d2bph_rr, d2bph_rt, d2bph_rp, d2bph_tt, d2bph_tp, d2bph_pp, G_test)
    
    ! Just check that we get reasonable values (not NaN or infinity)
    if (.not. (abs(sqg_test) < 1.0d20)) then
      print *, "ERROR: sqg_test is not finite:", sqg_test
      error stop 1
    end if
    
    if (.not. (abs(B_vartheta_test) < 1.0d20)) then
      print *, "ERROR: B_vartheta_test is not finite:", B_vartheta_test
      error stop 1
    end if
    
    if (.not. (abs(B_varphi_test) < 1.0d20)) then
      print *, "ERROR: B_varphi_test is not finite:", B_varphi_test
      error stop 1
    end if
    
    print *, "test_canonical_full_output: PASSED"
    
    ! Clean up
    call deallocate_can_coord()
    if (allocated(sA_phi)) deallocate(sA_phi)
    
  end subroutine test_canonical_full_output

end program test_canonical_full_regression