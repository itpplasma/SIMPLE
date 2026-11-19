program test_vmec_netcdf4
    !> A netCDF-4 (HDF5) copy of a 64-bit offset wout written by netCDF4-python
    !> must load into the same VMEC arrays. Its root group holds more links than
    !> one B-tree leaf, so this also covers the dense-link index of the reader.
    use new_vmec_stuff_mod, only: netcdffile, nstrm, kpar, rmnc, zmns, almns, &
                                  rmns, zmnc, almnc, aiota, phi, sps, axm, axn, s
    use vmec_alloc_sub, only: new_allocate_vmec_stuff, new_deallocate_vmec_stuff
    use vmecin_sub, only: vmecin

    implicit none

    double precision, allocatable :: ref_rmnc(:, :), ref_zmns(:, :), ref_almns(:, :)
    double precision, allocatable :: ref_aiota(:), ref_phi(:), ref_axm(:), ref_axn(:)
    double precision :: flux, ref_flux
    character(len=1000) :: classic_file, netcdf4_file

    call get_command_argument(1, classic_file)
    call get_command_argument(2, netcdf4_file)

    netcdffile = classic_file
    call new_allocate_vmec_stuff
    call vmecin(rmnc, zmns, almns, rmns, zmnc, almnc, aiota, phi, sps, axm, axn, s, &
                kpar + 1, nstrm, kpar, ref_flux)
    ref_rmnc = rmnc
    ref_zmns = zmns
    ref_almns = almns
    ref_aiota = aiota
    ref_phi = phi
    ref_axm = axm
    ref_axn = axn
    call new_deallocate_vmec_stuff

    netcdffile = netcdf4_file
    call new_allocate_vmec_stuff
    call vmecin(rmnc, zmns, almns, rmns, zmnc, almnc, aiota, phi, sps, axm, axn, s, &
                kpar + 1, nstrm, kpar, flux)

    if (any(shape(rmnc) /= shape(ref_rmnc))) error stop "netCDF-4 wout shape differs"
    if (any(rmnc /= ref_rmnc) .or. any(zmns /= ref_zmns) .or. any(almns /= ref_almns)) &
        error stop "netCDF-4 wout Fourier harmonics differ"
    if (any(aiota /= ref_aiota) .or. any(phi /= ref_phi) .or. flux /= ref_flux) &
        error stop "netCDF-4 wout profiles differ"
    if (any(axm /= ref_axm) .or. any(axn /= ref_axn)) &
        error stop "netCDF-4 wout mode numbers differ"
    call new_deallocate_vmec_stuff
    print *, "netCDF-4 wout matches the 64-bit offset wout"
end program test_vmec_netcdf4
