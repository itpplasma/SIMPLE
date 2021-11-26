!
  module chamb_mod
    logical :: rnegflag=.false.
!$omp threadprivate(rnegflag)
  end module chamb_mod
!
  module parmot_mod
    double precision :: rmu,ro0,eeff
  end module parmot_mod
!
  module new_vmec_stuff_mod
    character*1000   :: netcdffile
    integer          :: nsurfm,nstrm,nper,kpar
    integer          :: multharm,n_theta,n_phi
    integer          :: ns_A=5  !<- spline order for vector potential
    integer          :: ns_s=5  !<- spline order for R,Z,lambda over s
    integer          :: ns_tp=5 !<- spline order for R,Z,lambda over theta and varphi
    double precision :: rmajor,h_theta,h_phi
    double precision, dimension(:),     allocatable :: axm,axn,soa
    double precision, dimension(:),     allocatable :: aiota,s,sps,phi
    double precision, dimension(:,:),   allocatable :: almnc,rmnc,zmnc
    double precision, dimension(:,:),   allocatable :: almns,rmns,zmns
!
    double precision, dimension(:,:,:,:,:,:), allocatable :: sR,sZ,slam
  end module new_vmec_stuff_mod
!
  module vector_potentail_mod
    integer :: ns
    double precision :: hs,torflux
    double precision, dimension(:,:),         allocatable :: sA_phi
  end module vector_potentail_mod
!
  module canonical_coordinates_mod
    integer, parameter :: ns_max=6, n_qua=3
!
    integer :: ns_s_c,ns_tp_c
    integer :: ns_c,n_theta_c,n_phi_c,nh_stencil
    double precision :: hs_c,h_theta_c,h_phi_c
    double precision, dimension(:,:,:), allocatable :: G_c,sqg_c,B_vartheta_c,B_varphi_c
!
    double precision, dimension(ns_max)                     :: derf1,derf2,derf3
    double precision, dimension(:,:,:,:,:,:),   allocatable :: s_G_c
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: s_sqg_Bt_Bp
  end module canonical_coordinates_mod
!
  module boozer_coordinates_mod
    integer, parameter :: ns_max=6, n_qua=3
!
    logical :: use_B_r=.false., use_del_tp_B=.false.
!
    integer :: ns_s_B,ns_tp_B
    integer :: ns_B,n_theta_B,n_phi_B
    double precision :: hs_B,h_theta_B,h_phi_B
!
    double precision, dimension(ns_max)                     :: derf1,derf2,derf3
!
! spline coefficients for Boozer $B_\vartheta$ and $B_\varphi$:
    double precision, dimension(:,:,:),         allocatable :: s_Bcovar_tp_B
!
! spline coefficients for Boozer $B$ and $B_r$:
    double precision, dimension(:,:,:,:,:,:),   allocatable :: s_Bmod_B,s_Bcovar_r_B
!
! spline coefficients for $\Delta \vartheta_{BV}=\vartheta_B-\theta_V$
! and $\Delta \varphi_{BV}=\varphi_B-\varphi_V=G$ as functions of VMEC coordinates, s_delt_delp_V,
! and as functions of Boozer coordinates, s_delt_delp_B:
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: s_delt_delp_V,s_delt_delp_B
  end module boozer_coordinates_mod
!
  module velo_mod
    integer :: isw_field_type = 1
  end module velo_mod
!
module gbpi_mod
  character*24 :: filed
  integer :: ierrfield
end module gbpi_mod
!
module diag_mod
  logical :: dodiag=.false.
  integer(8) :: icounter
end module diag_mod
