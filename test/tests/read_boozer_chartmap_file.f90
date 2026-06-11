program read_boozer_chartmap_file
    use boozer_chartmap_io, only: boozer_chartmap_data_t, read_boozer_chartmap
    implicit none

    character(len=512) :: filename
    type(boozer_chartmap_data_t) :: data

    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) then
        print *, 'usage: read_boozer_chartmap_file <chartmap.nc>'
        error stop 2
    end if

    call read_boozer_chartmap(trim(filename), data)
    print *, 'read_boozer_chartmap_file: read succeeded'
end program read_boozer_chartmap_file
