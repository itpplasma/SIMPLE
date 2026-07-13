program test_cli_fast_class
    use params, only: class_plot, fast_class, ntcut
    use simple_main, only: classification_enabled
    implicit none

    ntcut = -1
    class_plot = .false.
    fast_class = .false.
    if (classification_enabled()) error stop "classification enabled without a mode"

    fast_class = .true.
    if (.not. classification_enabled()) error stop "fast classification not enabled"

    fast_class = .false.
    ntcut = 1
    if (.not. classification_enabled()) error stop "cut classification not enabled"

    ntcut = -1
    class_plot = .true.
    if (.not. classification_enabled()) error stop "plot classification not enabled"
end program test_cli_fast_class
