module nctools_module

  use netcdf
  implicit none

  interface nc_inq_dim
     module procedure nc_inq_dim_1
     module procedure nc_inq_dim_2
  end interface nc_inq_dim
  
  interface nc_get
     module procedure nc_get_int_0
     module procedure nc_get_int_1
     module procedure nc_get_double_0
     module procedure nc_get_double_1
     module procedure nc_get_double_2
     module procedure nc_get_double_3
  end interface nc_get
  
contains

  subroutine nc_get_int_0(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    integer, intent(out) :: var
    integer :: varid
    
    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_int_0

  subroutine nc_get_int_1(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    integer, dimension(:), intent(out) :: var
    integer :: varid
    
    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_int_1
  
  subroutine nc_get_double_0(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, intent(out) :: var
    integer :: varid
    
    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_0

  subroutine nc_get_double_1(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:), intent(out) :: var
    integer :: varid
    
    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_1

  subroutine nc_get_double_2(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:,:), intent(out) :: var
    integer :: varid
    
    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_2

  subroutine nc_get_double_3(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:,:,:), intent(out) :: var
    integer :: varid
    
    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_3
  
  subroutine nc_inq_dim_1(ncid, name, len)
    integer :: ncid, varid
    character(len=*)      :: name
    integer, dimension(1) :: dimids
    integer :: len

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_inquire_variable(ncid, varid, dimids = dimids))
    call nf90_check(nf90_inquire_dimension(ncid, dimids(1), len = len))
  end subroutine nc_inq_dim_1

  subroutine nc_inq_dim_2(ncid, name, len)
    integer :: ncid, varid
    character(len=*)      :: name
    integer, dimension(2) :: dimids
    integer, dimension(2) :: len

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_inquire_variable(ncid, varid, dimids = dimids))
    call nf90_check(nf90_inquire_dimension(ncid, dimids(1), len = len(1)))
    call nf90_check(nf90_inquire_dimension(ncid, dimids(2), len = len(2)))
  end subroutine nc_inq_dim_2
  
  subroutine nf90_check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       call abort
       stop
    end if
  end subroutine nf90_check

  subroutine nc_close(ncid)
    integer :: ncid
    
    call nf90_check(nf90_close(ncid))
  end subroutine nc_close

  subroutine nc_open(filename, ncid)
    character(len=*) :: filename
    integer :: ncid

    call nf90_check(nf90_open(filename, NF90_NOWRITE, ncid))
  end subroutine nc_open

end module nctools_module
