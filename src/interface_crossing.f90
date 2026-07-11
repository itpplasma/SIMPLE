module interface_crossing
    !> Level-0 guiding-center interface crossing for SPECTRE stacked-rho charts
    !> (#443). apply_crossing is the single entry point; the Level-1 refraction
    !> map (#440) later replaces the map behind it without touching callers, so
    !> the level argument selects the map (only level 0 is implemented here).
    !>
    !> Level 0 holds (theta, zeta, mu), keeps the sign of v_par, switches volume,
    !> and rescales v_par from exact energy conservation across the pressure jump
    !> [[B]]. It is energy-exact in both branches (crossing and reflection).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: magfie

    implicit none
    private

    public :: apply_crossing, crossing_info_t, axis_offset
    public :: crossing_log_reset, crossing_log_record, crossing_log_write, &
              crossing_log_count
    public :: CROSSING_LEVEL0, CROSS_CROSSING, CROSS_REFLECTION, CROSS_LOSS, &
              CROSS_STOP

    integer, parameter :: CROSSING_LEVEL0 = 0
    integer, parameter :: CROSS_CROSSING = 1
    integer, parameter :: CROSS_REFLECTION = 2
    integer, parameter :: CROSS_LOSS = 3
    !> Symplectic orbits terminate at volume boundaries until the exact-landing
    !> crossing is wired (SIMPLE#441); their stop events share this log.
    integer, parameter :: CROSS_STOP = 4

    !> Inner cutoff for the innermost volume: rho_g = 0 is the coordinate axis
    !> where sqrt(g) = 0, so the marker reflects trivially before reaching it
    !> rather than crossing a (non-existent) inner interface.
    real(dp), parameter :: axis_offset = 1.0d-3

    !> Both-side |B| is read at rho_g = k -+ iface_eps: int(k -+ 1e-12) selects
    !> the neighbouring volume k-1 resp. k, and the libneo polynomial basis of
    !> each volume is well-defined at that offset from the shared interface point.
    real(dp), parameter :: iface_eps = 1.0d-12
    !> After the map rho_g is nudged reflect_nudge off the interface so the
    !> restarted RK45 step does not immediately re-detect the same event. It must
    !> clear the odeint event tolerance (~sqrt(epsilon) ~ 1.5e-8) or the next
    !> substep triggers a zero-length event at its start point.
    real(dp), parameter :: reflect_nudge = 1.0d-6

    type :: crossing_info_t
        integer :: event_type = 0
        integer :: iface = 0
        integer :: vol_from = 0
        integer :: vol_to = 0
        real(dp) :: theta = 0.0_dp
        real(dp) :: zeta = 0.0_dp
        real(dp) :: vpar_before = 0.0_dp
        real(dp) :: vpar_after = 0.0_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: bmod_home = 0.0_dp
        real(dp) :: bmod_target = 0.0_dp
    end type crossing_info_t

    type :: crossing_record_t
        integer :: ipart = 0
        real(dp) :: time = 0.0_dp
        type(crossing_info_t) :: info
    end type crossing_record_t

    type(crossing_record_t), allocatable :: crossing_log(:)
    integer :: crossing_count = 0

contains

    subroutine apply_crossing(y_iface, iface, direction, mvol, level, y_out, info)
        !> Map the guiding-center state y = (rho_g, theta, zeta, p, lambda) landed
        !> on interface rho_g = iface across to the neighbour volume. direction is
        !> +1 outward (home volume below the interface) or -1 inward (home above).
        real(dp), intent(in) :: y_iface(5)
        integer, intent(in) :: iface, direction, mvol, level
        real(dp), intent(out) :: y_out(5)
        type(crossing_info_t), intent(out) :: info

        real(dp) :: rho_face, rho_home, rho_target
        real(dp) :: bmod_home, bmod_target, perp_inv, mu, vpar, radicand

        if (level /= CROSSING_LEVEL0) then
            error stop 'apply_crossing: Level-1 crossing not yet implemented'
        end if

        rho_face = real(iface, dp)
        y_out = y_iface
        vpar = y_iface(4)*y_iface(5)

        info%iface = iface
        info%theta = y_iface(2)
        info%zeta = y_iface(3)
        info%vpar_before = vpar
        info%vpar_after = vpar

        if (direction == 1) then
            info%vol_from = iface
            info%vol_to = iface + 1
            rho_home = rho_face - iface_eps
            rho_target = rho_face + iface_eps
        else
            info%vol_from = iface + 1
            info%vol_to = iface
            rho_home = rho_face + iface_eps
            rho_target = rho_face - iface_eps
        end if

        bmod_home = bmod_at([rho_home, y_iface(2), y_iface(3)])
        info%bmod_home = bmod_home
        info%bmod_target = bmod_home

        ! perp_inv = v_perp^2/|B| = z(4)^2 (1 - z(5)^2)/|B| is the perpendicular
        ! invariant the RK45 path stores; mu = perp_inv/2 is the magnetic moment
        ! held fixed across the interface, so the whole map runs on this invariant.
        perp_inv = y_iface(4)**2*(1.0_dp - y_iface(5)**2)/bmod_home
        mu = 0.5_dp*perp_inv
        info%mu = mu

        if (direction == 1 .and. iface == mvol) then
            info%event_type = CROSS_LOSS
            info%vol_to = mvol + 1
            return
        end if

        if (direction == -1 .and. iface == 0) then
            info%event_type = CROSS_REFLECTION
            info%vol_to = info%vol_from
            info%vpar_after = -vpar
            y_out(5) = -y_iface(5)
            y_out(1) = axis_offset + reflect_nudge
            return
        end if

        bmod_target = bmod_at([rho_target, y_iface(2), y_iface(3)])
        info%bmod_target = bmod_target

        radicand = vpar**2 - 2.0_dp*mu*(bmod_target - bmod_home)
        if (radicand >= 0.0_dp) then
            info%event_type = CROSS_CROSSING
            info%vpar_after = sign(sqrt(radicand), vpar)
            y_out(5) = info%vpar_after/y_iface(4)
            y_out(1) = rho_face + real(direction, dp)*reflect_nudge
        else
            ! Forbidden crossing into the higher-field volume is a magnetic
            ! mirror: the marker stays home with v_par reversed (the same-side
            ! energy root), so |v_par|, mu and |B| are unchanged and H is exactly
            ! conserved while the orbit streams back off the sheet.
            info%event_type = CROSS_REFLECTION
            info%vol_to = info%vol_from
            info%vpar_after = -vpar
            y_out(5) = -y_iface(5)
            y_out(1) = rho_face - real(direction, dp)*reflect_nudge
        end if
    end subroutine apply_crossing

    function bmod_at(x) result(bmod)
        real(dp), intent(in) :: x(3)
        real(dp) :: bmod

        real(dp) :: sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end function bmod_at

    subroutine crossing_log_reset(capacity)
        integer, intent(in) :: capacity

        if (allocated(crossing_log)) deallocate (crossing_log)
        allocate (crossing_log(max(capacity, 1)))
        crossing_count = 0
    end subroutine crossing_log_reset

    subroutine crossing_log_record(ipart, time, info)
        integer, intent(in) :: ipart
        real(dp), intent(in) :: time
        type(crossing_info_t), intent(in) :: info

        type(crossing_record_t), allocatable :: tmp(:)

        !$omp critical (spectre_crossing_log)
        if (.not. allocated(crossing_log)) allocate (crossing_log(64))
        if (crossing_count >= size(crossing_log)) then
            allocate (tmp(2*size(crossing_log)))
            tmp(1:crossing_count) = crossing_log(1:crossing_count)
            call move_alloc(tmp, crossing_log)
        end if
        crossing_count = crossing_count + 1
        crossing_log(crossing_count)%ipart = ipart
        crossing_log(crossing_count)%time = time
        crossing_log(crossing_count)%info = info
        !$omp end critical (spectre_crossing_log)
    end subroutine crossing_log_record

    integer function crossing_log_count()
        crossing_log_count = crossing_count
    end function crossing_log_count

    subroutine crossing_log_write(filename)
        character(*), intent(in) :: filename

        integer :: unit, i
        integer, allocatable :: order(:)

        if (.not. allocated(crossing_log)) return

        ! Checkpoint dumps call this from a worker thread while others append, so
        ! serialize against crossing_log_record to keep the growable array live.
        !$omp critical (spectre_crossing_log)
        call sorted_by_particle(order)
        open (newunit=unit, file=filename, recl=1024)
        write (unit, '(A)') '# particle  time  iface  type  vol_from  vol_to  '// &
            'theta  zeta  vpar_before  vpar_after  mu  bmod_home  bmod_target'
        do i = 1, crossing_count
            associate (r => crossing_log(order(i)))
                write (unit, *) r%ipart, r%time, r%info%iface, r%info%event_type, &
                    r%info%vol_from, r%info%vol_to, r%info%theta, r%info%zeta, &
                    r%info%vpar_before, r%info%vpar_after, r%info%mu, &
                    r%info%bmod_home, r%info%bmod_target
            end associate
        end do
        close (unit)
        !$omp end critical (spectre_crossing_log)
    end subroutine crossing_log_write

    subroutine sorted_by_particle(order)
        !> Stable counting sort of the log by particle index. Each particle is
        !> traced by a single thread, so its events are appended in time order;
        !> grouping by particle therefore yields a thread-order-independent,
        !> reproducible file whose per-particle blocks stay time-ordered.
        integer, allocatable, intent(out) :: order(:)

        integer :: i, p, pmax
        integer, allocatable :: cnt(:), pos(:)

        allocate (order(crossing_count))
        if (crossing_count == 0) return

        pmax = 0
        do i = 1, crossing_count
            pmax = max(pmax, crossing_log(i)%ipart)
        end do

        allocate (cnt(pmax), pos(pmax))
        cnt = 0
        do i = 1, crossing_count
            p = crossing_log(i)%ipart
            cnt(p) = cnt(p) + 1
        end do
        pos(1) = 1
        do p = 2, pmax
            pos(p) = pos(p - 1) + cnt(p - 1)
        end do
        do i = 1, crossing_count
            p = crossing_log(i)%ipart
            order(pos(p)) = i
            pos(p) = pos(p) + 1
        end do
    end subroutine sorted_by_particle

end module interface_crossing
