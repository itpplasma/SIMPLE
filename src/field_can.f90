module field_can_mod
use diag_mod, only : icounter
use boozer_sub, only : splint_boozer_coord
use magfie_sub, only : TEST, CANFLUX, BOOZER, MEISS, ALBERT
use field, only : magnetic_field_t, vmec_field_t
use field_can_base, only : twopi, evaluate_base => evaluate, coordinate_transform, &
  identity_transform, FieldCan
use field_can_test, only : evaluate_test
use field_can_flux, only : evaluate_flux, integ_to_ref_flux, ref_to_integ_flux
use field_can_boozer, only : evaluate_boozer, integ_to_ref_boozer, ref_to_integ_boozer
use field_can_meiss, only : init_meiss, evaluate_meiss, &
  integ_to_ref_meiss, ref_to_integ_meiss
use field_can_albert, only : evaluate_albert, init_albert, integ_to_ref_albert, &
  ref_to_integ_albert

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)

procedure(evaluate_base), pointer :: evaluate => null()

  ! Conversion to and from reference coordinates - currently VMEC coordinates (s, th, ph)
procedure(coordinate_transform), pointer :: integ_to_ref => identity_transform
procedure(coordinate_transform), pointer :: ref_to_integ => identity_transform

contains

subroutine field_can_from_name(field_name, field_noncan)

  character(*), intent(in) :: field_name
  !> For FieldCanMeiss
  class(magnetic_field_t), intent(in), optional :: field_noncan

  select case(trim(field_name))
    case("test")
      evaluate => evaluate_test
    case("flux")
      evaluate => evaluate_flux
      integ_to_ref => integ_to_ref_flux
      ref_to_integ => ref_to_integ_flux
    case("boozer")
      evaluate => evaluate_boozer
      integ_to_ref => integ_to_ref_boozer
      ref_to_integ => ref_to_integ_boozer
    case("meiss")
      if (present(field_noncan)) then
        call init_meiss(field_noncan)
      end if
      evaluate => evaluate_meiss
      integ_to_ref => integ_to_ref_meiss
      ref_to_integ => ref_to_integ_meiss
    case("albert")
      if (present(field_noncan)) then
        call init_albert(field_noncan)
      end if
      evaluate => evaluate_albert
      integ_to_ref => integ_to_ref_albert
      ref_to_integ => ref_to_integ_albert
    case default
      print *, "field_can_from_name: Unknown field type ", field_name
      error stop
  end select
end subroutine field_can_from_name


subroutine field_can_from_id(field_id, field_noncan)
  integer, intent(in) :: field_id
  !> For FieldCanMeiss
  class(magnetic_field_t), intent(in), optional :: field_noncan

  if (present(field_noncan)) then
    call field_can_from_name(name_from_id(field_id), field_noncan)
  else
    call field_can_from_name(name_from_id(field_id))
  end if

end subroutine field_can_from_id


function name_from_id(field_id)
  character(128) :: name_from_id
  integer, intent(in) :: field_id

  select case(field_id)
    case(TEST)
      name_from_id = "test"
    case(CANFLUX)
      name_from_id = "flux"
    case(BOOZER)
      name_from_id = "boozer"
    case(MEISS)
      name_from_id = "meiss"
    case(ALBERT)
      name_from_id = "albert"
    case default
      print *, "name_from_id: Unknown field id ", field_id
      error stop
  end select
end function name_from_id


function id_from_name(field_name)
  integer :: id_from_name
  character(*), intent(in) :: field_name

  select case(trim(field_name))
    case("test")
      id_from_name = TEST
    case("flux")
      id_from_name = CANFLUX
    case("boozer")
      id_from_name = BOOZER
    case("meiss")
      id_from_name = MEISS
    case("albert")
      id_from_name = ALBERT
    case default
      print *, "id_from_name: Unknown field type ", field_name
      error stop
  end select
end function id_from_name


subroutine init_field_can(field_id, field_noncan)
  use get_can_sub, only : get_canonical_coordinates, get_canonical_coordinates_with_field
  use boozer_sub, only : get_boozer_coordinates, get_boozer_coordinates_with_field
  use field_can_meiss, only : get_meiss_coordinates
  use field_can_albert, only : get_albert_coordinates

  integer, intent(in) :: field_id
  class(magnetic_field_t), intent(in), optional :: field_noncan
  class(magnetic_field_t), allocatable :: field_to_use

  if (present(field_noncan)) then
    allocate(field_to_use, source=field_noncan)
    call field_can_from_id(field_id, field_noncan)
  else
    allocate(field_to_use, source=vmec_field_t())
    call field_can_from_id(field_id, vmec_field_t())
  end if
  
  select case (field_id)
    case (TEST)
      call field_can_from_id(field_id, field_to_use)
    case (CANFLUX)
      call get_canonical_coordinates_with_field(field_to_use)
    case (BOOZER)
      call get_boozer_coordinates_with_field(field_to_use)
    case (MEISS)
      call get_meiss_coordinates
    case (ALBERT)
      call get_albert_coordinates
    case default
      print *, "init_field_can: Unknown field id ", field_id
      error stop
  end select
  
  deallocate(field_to_use)
end subroutine init_field_can


subroutine FieldCan_init(f, mu, ro0, vpar)
  type(FieldCan), intent(inout) :: f
  real(dp), intent(in), optional  :: mu, ro0, vpar

  if (present(mu)) then
    f%mu = mu
  else
    f%mu = 0d0
  end if

  if (present(ro0)) then
    f%ro0 = ro0
  else
    f%ro0 = 0d0
  end if

  if (present(vpar)) then
    f%vpar = vpar
  else
    f%vpar = 0d0
  end if

end subroutine FieldCan_init


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine get_val(f, pphi)
  !
  ! computes values of H, pth and vpar at z=(r, th, ph, pphi)
  !
  !
  type(FieldCan), intent(inout) :: f
  real(dp), intent(in) :: pphi

  f%vpar = (pphi - f%Aph/f%ro0)/f%hph
  f%H = f%vpar**2/2d0 + f%mu*f%Bmod
  f%pth = f%hth*f%vpar + f%Ath/f%ro0

end subroutine get_val


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine get_derivatives(f, pphi)
  !
  ! computes H, pth and vpar at z=(r, th, ph, pphi) and their derivatives
  !
  !
  type(FieldCan), intent(inout) :: f
  real(dp), intent(in) :: pphi

  call get_val(f, pphi)

  f%dvpar(1:3) = -(f%dAph/f%ro0 + f%dhph*f%vpar)/f%hph
  f%dvpar(4)   = 1d0/f%hph

  f%dH(1:3) = f%vpar*f%dvpar(1:3) + f%mu*f%dBmod
  f%dH(4)   = f%vpar/f%hph

  f%dpth(1:3) = f%dvpar(1:3)*f%hth + f%vpar*f%dhth + f%dAth/f%ro0

  f%dpth(4) = f%hth/f%hph

end subroutine get_derivatives

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine get_derivatives2(f, pphi)
  !
  ! computes H, pth and vpar at z=(r, th, ph, pphi) up to 2nd derivatives
  ! order of second derivatives:
  ! d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
  ! d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
  !
  type(FieldCan), intent(inout) :: f
  real(dp), intent(in) :: pphi

  call get_derivatives(f, pphi)

  f%d2vpar(1:6) = -f%d2Aph/f%ro0 - f%d2hph*f%vpar
  f%d2vpar(1) = f%d2vpar(1) - 2d0*f%dhph(1)*f%dvpar(1)
  f%d2vpar(2) = f%d2vpar(2) - (f%dhph(1)*f%dvpar(2) + f%dhph(2)*f%dvpar(1))
  f%d2vpar(3) = f%d2vpar(3) - (f%dhph(1)*f%dvpar(3) + f%dhph(3)*f%dvpar(1))
  f%d2vpar(4) = f%d2vpar(4) - 2d0*f%dhph(2)*f%dvpar(2)
  f%d2vpar(5) = f%d2vpar(5) - (f%dhph(2)*f%dvpar(3) + f%dhph(3)*f%dvpar(2))
  f%d2vpar(6) = f%d2vpar(6) - 2d0*f%dhph(3)*f%dvpar(3)
  f%d2vpar(1:6) = f%d2vpar(1:6)/f%hph

  f%d2H(1:6) = f%vpar*f%d2vpar(1:6) + f%mu*f%d2Bmod ! + qi*d2Phie
  f%d2H(1) = f%d2H(1) + f%dvpar(1)**2
  f%d2H(2) = f%d2H(2) + f%dvpar(1)*f%dvpar(2)
  f%d2H(3) = f%d2H(3) + f%dvpar(1)*f%dvpar(3)
  f%d2H(4) = f%d2H(4) + f%dvpar(2)**2
  f%d2H(5) = f%d2H(5) + f%dvpar(2)*f%dvpar(3)
  f%d2H(6) = f%d2H(6) + f%dvpar(3)**2

  f%d2pth(1:6) = f%d2vpar(1:6)*f%hth + f%vpar*f%d2hth + f%d2Ath/f%ro0
  f%d2pth(1) = f%d2pth(1) + 2d0*f%dvpar(1)*f%dhth(1)
  f%d2pth(2) = f%d2pth(2) + f%dvpar(1)*f%dhth(2) + f%dvpar(2)*f%dhth(1)
  f%d2pth(3) = f%d2pth(3) + f%dvpar(1)*f%dhth(3) + f%dvpar(3)*f%dhth(1)
  f%d2pth(4) = f%d2pth(4) + 2d0*f%dvpar(2)*f%dhth(2)
  f%d2pth(5) = f%d2pth(5) + f%dvpar(2)*f%dhth(3) + f%dvpar(3)*f%dhth(2)
  f%d2pth(6) = f%d2pth(6) + 2d0*f%dvpar(3)*f%dhth(3)

  f%d2vpar(7:9) = -f%dhph/f%hph**2
  f%d2H(7:9) = f%dvpar(1:3)/f%hph + f%vpar*f%d2vpar(7:9)
  f%d2pth(7:9) = f%dhth/f%hph + f%hth*f%d2vpar(7:9)

end subroutine get_derivatives2

! Wrapper subroutines to expose function pointers to Python
subroutine integ_to_ref_wrapper(xinteg, xref)
  real(dp), intent(in) :: xinteg(3)
  real(dp), intent(out) :: xref(3)
  call integ_to_ref(xinteg, xref)
end subroutine integ_to_ref_wrapper

subroutine ref_to_integ_wrapper(xref, xinteg)
  real(dp), intent(in) :: xref(3)
  real(dp), intent(out) :: xinteg(3)
  call ref_to_integ(xref, xinteg)
end subroutine ref_to_integ_wrapper

end module field_can_mod
