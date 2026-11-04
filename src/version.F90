module version
    !> SIMPLE version information. CMake defines SIMPLE_VERSION_STRING from
    !> git describe; builds without it, such as fpm, report "unknown".
    !>
    !> Version format: <tag>-<commits>-g<hash>[-dirty]
    !> Example: v1.5.2-42-g62fed6e or v1.5.2-42-g62fed6e-dirty
    implicit none

#ifdef SIMPLE_VERSION_STRING
    character(len=*), parameter :: simple_version = SIMPLE_VERSION_STRING
#else
    character(len=*), parameter :: simple_version = "unknown"
#endif

end module version
