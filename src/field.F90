module field
    !> Field module aggregating all field types and factory functions.

    use libneo_coordinates, only: detect_refcoords_file_type, refcoords_file_chartmap, &
                                  refcoords_file_vmec_wout, refcoords_file_unknown
    use field_base, only: magnetic_field_t
    use field_vmec, only: vmec_field_t, create_vmec_field
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use field_boozer_chartmap, only: boozer_chartmap_field_t, &
                                     create_boozer_chartmap_field, &
                                     is_boozer_chartmap
    use field_spectre, only: spectre_field_t, create_spectre_field

    implicit none

contains

    subroutine field_clone(source, dest)
        class(magnetic_field_t), intent(in) :: source
        class(magnetic_field_t), allocatable, intent(out) :: dest

        select type (source)
        type is (vmec_field_t)
            allocate (vmec_field_t :: dest)
            select type (dest)
            type is (vmec_field_t)
                dest = source
            class default
                error stop 'field_clone: Allocation failure (vmec)'
            end select
        type is (coils_field_t)
            allocate (coils_field_t :: dest)
            select type (dest)
            type is (coils_field_t)
                dest = source
            class default
                error stop 'field_clone: Allocation failure (coils)'
            end select
        type is (splined_field_t)
            allocate (splined_field_t :: dest)
            select type (dest)
            type is (splined_field_t)
                dest = source
            class default
                error stop 'field_clone: Allocation failure (splined)'
            end select
        type is (boozer_chartmap_field_t)
            allocate (boozer_chartmap_field_t :: dest)
            select type (dest)
            type is (boozer_chartmap_field_t)
                dest = source
            class default
                error stop 'field_clone: Allocation failure (boozer_chartmap)'
            end select
        type is (spectre_field_t)
            allocate (spectre_field_t :: dest)
            select type (dest)
            type is (spectre_field_t)
                dest = source
            class default
                error stop 'field_clone: Allocation failure (spectre)'
            end select
        class default
            error stop 'field_clone: Unsupported field type'
        end select
    end subroutine field_clone

    subroutine field_from_file(filename, field)
        !> Create the field from a file, chosen by the file's content
        !> (field_file_kind). Coil files and fields on a chartmap become a
        !> splined_field_t on the reference coordinates.
        use reference_coordinates, only: ref_coords
        ! Imported here, not by the module: fortio's IEEE symbols would make
        ! gfortran save and restore the FP environment in every procedure of
        ! each module that uses field.
        use field_chartmap_nodes, only: create_splined_field_from_chartmap
        use libneo_chartmap_field, only: chartmap_has_field
        use field_file_kind, only: field_file_kind_of, hdf5_root_members, &
                                   FIELD_FILE_FROZEN_CURRENT, FIELD_FILE_SPECTRE, &
                                   FIELD_FILE_HDF5, FIELD_FILE_NETCDF, FIELD_FILE_COILS

        character(*), intent(in) :: filename
        class(magnetic_field_t), allocatable, intent(out) :: field

        type(coils_field_t) :: raw_coils
        type(splined_field_t), allocatable :: splined_coils
        type(vmec_field_t) :: vmec_field
        type(boozer_chartmap_field_t), allocatable :: bc_temp
        type(spectre_field_t), allocatable :: spectre_temp
        integer :: kind, file_type, ierr
        character(len=2048) :: message

        kind = field_file_kind_of(filename)
        select case (kind)
        case (FIELD_FILE_FROZEN_CURRENT)
            print *, 'field_from_file: ', trim(filename), ' is a solved ', &
                'frozen-current field; evaluate it on the chart once with libneo ', &
                'frozen_current_chartmap.x and use that file as field_input'
            error stop
        case (FIELD_FILE_SPECTRE)
            allocate (spectre_temp)
            call create_spectre_field(spectre_temp, filename, ierr)
            if (ierr /= 0) then
                print *, 'field_from_file: create_spectre_field failed for ', &
                    trim(filename), ' (ierr = ', ierr, ')'
                error stop
            end if
            call move_alloc(spectre_temp, field)
        case (FIELD_FILE_HDF5, FIELD_FILE_NETCDF)
            if (chartmap_has_field(filename)) then
                allocate (splined_coils)
                call create_splined_field_from_chartmap(filename, ref_coords, &
                                                        splined_coils)
                call move_alloc(splined_coils, field)
                return
            end if
            if (is_boozer_chartmap(filename)) then
                call create_boozer_chartmap_field(filename, bc_temp)
                call move_alloc(bc_temp, field)
                return
            end if
            call detect_refcoords_file_type(filename, file_type, ierr, message)
            if (ierr /= 0) file_type = refcoords_file_unknown
            select case (file_type)
            case (refcoords_file_vmec_wout)
                call create_vmec_field(vmec_field)
                call field_clone(vmec_field, field)
            case (refcoords_file_chartmap)
                print *, &
                    'field_from_file: chartmap NetCDF is a coordinate system ', &
                    'file, not a field:'
                print *, '  filename = ', trim(filename)
                print *, &
                    'Set coord_input to this chartmap file and set field_input', &
                    ' to a VMEC wout or to a field on this chartmap.'
                error stop
            case default
                if (kind == FIELD_FILE_HDF5) then
                    print *, 'field_from_file: ', trim(filename), ' is HDF5 but ', &
                        'neither a field on a chartmap, a frozen-current field ', &
                        '(u, kmax, coil_points), a SPECTRE equilibrium ', &
                        '(input/physics/Nvol, output/Rbc), a VMEC wout nor a ', &
                        'Boozer chartmap'
                    print *, '  root members: ', hdf5_root_members(filename)
                else
                    print *, 'field_from_file: unknown NetCDF file type: ', &
                        trim(filename)
                    if (ierr /= 0) print *, trim(message)
                end if
                error stop
            end select
        case (FIELD_FILE_COILS)
            call create_coils_field(filename, raw_coils)
            allocate (splined_coils)
            call create_splined_field(raw_coils, ref_coords, splined_coils)
            call move_alloc(splined_coils, field)
        case default
            print *, 'field_from_file: ', trim(filename), ' is neither HDF5, ', &
                'NetCDF nor a coil file (point count, then x y z current rows)'
            error stop
        end select
    end subroutine field_from_file

    !> SPECTRE equilibrium by content: an HDF5 file with SPECTRE's input and
    !> output groups.
    logical function is_spectre_file(filename) result(is_spectre)
        use field_file_kind, only: field_file_kind_of, FIELD_FILE_SPECTRE
        character(*), intent(in) :: filename

        is_spectre = field_file_kind_of(filename) == FIELD_FILE_SPECTRE
    end function is_spectre_file

end module field
