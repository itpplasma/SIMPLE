program test_particle_rng
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use params, only: deterministic, ran_seed, reset_seed_if_deterministic, &
        particle_stream_id, particle_index_offset, startmode, idx, &
        ntestpart, batch_size, reuse_batch, init_batch
    use samplers, only: sample
    use omp_lib, only: omp_get_num_threads
    implicit none

    integer, parameter :: population = 16, draws = 12, subset = 5
    real(dp) :: serial(draws, population), parallel(draws, population)
    real(dp) :: replay(draws), legacy(draws), row(5), selected(5, subset)
    integer :: i, j, unit, seedsize, saved_idx(subset), observed_threads
    integer, allocatable :: seed(:)

    deterministic = .true.
    ran_seed = 314159
    call reset_seed_if_deterministic
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed)
    if (any(seed /= ran_seed)) error stop 'Legacy untagged seed changed'
    call random_number(legacy)
    call reset_seed_if_deterministic
    call random_number(replay)
    if (any(legacy /= replay)) error stop 'Legacy untagged replay changed'

    startmode = 2
    particle_index_offset = 100
    !$omp parallel do num_threads(1) private(i)
    do i = 1, population
        call reset_seed_if_deterministic(particle_stream_id(i))
        call random_number(serial(:, i))
    end do
    !$omp end parallel do
    observed_threads = 0
    !$omp parallel num_threads(4) private(i)
    !$omp single
    observed_threads = omp_get_num_threads()
    !$omp end single
    !$omp do schedule(dynamic, 1)
    do i = 1, population
        call reset_seed_if_deterministic(particle_stream_id(i))
        call random_number(parallel(:, i))
    end do
    !$omp end do
    !$omp end parallel
    if (observed_threads /= 4) error stop 'Four-thread RNG check not executed'
    if (any(serial /= parallel)) error stop 'Marker stream depends on thread count'
    do i = 1, population
        do j = i + 1, population
            if (all(serial(:, i) == serial(:, j))) &
                error stop 'Different markers share identical draws'
        end do
    end do

    ! A contiguous external shard retains the same original source-row IDs.
    particle_index_offset = 108
    do i = 1, population/2
        call reset_seed_if_deterministic(particle_stream_id(i))
        call random_number(replay)
        if (any(replay /= serial(:, i + 8))) &
            error stop 'External shard changed a global marker stream'
    end do

    ! Native selection persists row numbers as text and reads those exact rows.
    open(newunit=unit, file='start.dat', status='replace', action='write')
    do i = 1, population
        row = real(i, dp)
        write(unit, *) row
    end do
    close(unit)
    particle_index_offset = 100
    startmode = 4
    ntestpart = population
    batch_size = subset
    reuse_batch = .false.
    call init_batch
    if (ntestpart /= subset) error stop 'Native batch size not applied'
    if (any(idx < 1) .or. any(idx > population)) error stop 'Invalid source row'
    do i = 2, subset
        if (idx(i) <= idx(i - 1)) error stop 'Duplicate or unsorted source row'
    end do
    saved_idx = idx
    call sample(selected, reuse_batch)
    do i = 1, subset
        if (any(selected(:, i) /= real(idx(i), dp))) &
            error stop 'Native batch loaded the wrong source row'
        call reset_seed_if_deterministic(particle_stream_id(i))
        call random_number(replay)
        if (any(replay /= serial(:, idx(i)))) &
            error stop 'Native selection changed a global marker stream'
    end do

    ntestpart = population
    reuse_batch = .true.
    call init_batch
    if (any(idx /= saved_idx)) error stop 'Text batch indices did not round trip'
    call sample(selected, reuse_batch)
    do i = 1, subset
        if (any(selected(:, i) /= real(saved_idx(i), dp))) &
            error stop 'Reused batch loaded the wrong source row'
    end do
    ! Check both source-file endpoints against an independently specified list.
    saved_idx = [1, 3, 6, 12, population]
    open(newunit=unit, file='batch.dat', status='replace', action='write')
    do i = 1, subset
        write(unit, *) saved_idx(i)
    end do
    close(unit)
    ntestpart = population
    call init_batch
    call sample(selected, reuse_batch)
    do i = 1, subset
        if (any(selected(:, i) /= real(saved_idx(i), dp))) &
            error stop 'Explicit one-based batch endpoints loaded incorrectly'
        call reset_seed_if_deterministic(particle_stream_id(i))
        call random_number(replay)
        if (any(replay /= serial(:, saved_idx(i)))) &
            error stop 'Explicit batch changed original marker identity'
    end do
    print *, 'PASS: legacy replay, marker streams, threads, shards and batches'
end program test_particle_rng
