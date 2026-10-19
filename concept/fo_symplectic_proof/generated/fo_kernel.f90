module fo_kernel
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none
contains
subroutine fo_eval(p,aj,gj,mass,charge,e0,tj,energy,kenergy,gh,hh,f,df)
real(dp),intent(in)::p(3),aj(3,10),gj(6,10),mass,charge,e0,tj(10)
real(dp),intent(out)::energy,kenergy,gh(6),hh(6,6),f(6),df(6,6)
call fo_scalar( &
p(1), &
p(2), &
p(3), &
mass, charge, e0, &
aj(1,1), &
aj(1,2), &
aj(1,3), &
aj(1,4), &
aj(1,5), &
aj(1,6), &
aj(1,7), &
aj(1,8), &
aj(1,9), &
aj(1,10), &
aj(2,1), &
aj(2,2), &
aj(2,3), &
aj(2,4), &
aj(2,5), &
aj(2,6), &
aj(2,7), &
aj(2,8), &
aj(2,9), &
aj(2,10), &
aj(3,1), &
aj(3,2), &
aj(3,3), &
aj(3,4), &
aj(3,5), &
aj(3,6), &
aj(3,7), &
aj(3,8), &
aj(3,9), &
aj(3,10), &
gj(1,1), &
gj(1,2), &
gj(1,3), &
gj(1,4), &
gj(1,5), &
gj(1,6), &
gj(1,7), &
gj(1,8), &
gj(1,9), &
gj(1,10), &
gj(2,1), &
gj(2,2), &
gj(2,3), &
gj(2,4), &
gj(2,5), &
gj(2,6), &
gj(2,7), &
gj(2,8), &
gj(2,9), &
gj(2,10), &
gj(3,1), &
gj(3,2), &
gj(3,3), &
gj(3,4), &
gj(3,5), &
gj(3,6), &
gj(3,7), &
gj(3,8), &
gj(3,9), &
gj(3,10), &
gj(4,1), &
gj(4,2), &
gj(4,3), &
gj(4,4), &
gj(4,5), &
gj(4,6), &
gj(4,7), &
gj(4,8), &
gj(4,9), &
gj(4,10), &
gj(5,1), &
gj(5,2), &
gj(5,3), &
gj(5,4), &
gj(5,5), &
gj(5,6), &
gj(5,7), &
gj(5,8), &
gj(5,9), &
gj(5,10), &
gj(6,1), &
gj(6,2), &
gj(6,3), &
gj(6,4), &
gj(6,5), &
gj(6,6), &
gj(6,7), &
gj(6,8), &
gj(6,9), &
gj(6,10), &
tj(1), &
tj(2), &
tj(3), &
tj(4), &
tj(5), &
tj(6), &
tj(7), &
tj(8), &
tj(9), &
tj(10), &
energy,kenergy, &
gh(1), &
gh(2), &
gh(3), &
gh(4), &
gh(5), &
gh(6), &
hh(1,1), &
hh(2,1), &
hh(3,1), &
hh(4,1), &
hh(5,1), &
hh(6,1), &
hh(1,2), &
hh(2,2), &
hh(3,2), &
hh(4,2), &
hh(5,2), &
hh(6,2), &
hh(1,3), &
hh(2,3), &
hh(3,3), &
hh(4,3), &
hh(5,3), &
hh(6,3), &
hh(1,4), &
hh(2,4), &
hh(3,4), &
hh(4,4), &
hh(5,4), &
hh(6,4), &
hh(1,5), &
hh(2,5), &
hh(3,5), &
hh(4,5), &
hh(5,5), &
hh(6,5), &
hh(1,6), &
hh(2,6), &
hh(3,6), &
hh(4,6), &
hh(5,6), &
hh(6,6), &
f(1), &
f(2), &
f(3), &
f(4), &
f(5), &
f(6), &
df(1,1), &
df(2,1), &
df(3,1), &
df(4,1), &
df(5,1), &
df(6,1), &
df(1,2), &
df(2,2), &
df(3,2), &
df(4,2), &
df(5,2), &
df(6,2), &
df(1,3), &
df(2,3), &
df(3,3), &
df(4,3), &
df(5,3), &
df(6,3), &
df(1,4), &
df(2,4), &
df(3,4), &
df(4,4), &
df(5,4), &
df(6,4), &
df(1,5), &
df(2,5), &
df(3,5), &
df(4,5), &
df(5,5), &
df(6,5), &
df(1,6), &
df(2,6), &
df(3,6), &
df(4,6), &
df(5,6), &
df(6,6))
end subroutine fo_eval
    subroutine fo_eval_first(p,aj,gj,mass,charge,e0,tj,energy,kenergy,gh,hh,f,df)
        real(dp), intent(in)::p(3), aj(3, 10), gj(6, 10), mass, charge, e0, tj(10)
        real(dp), intent(out)::energy, kenergy, gh(6), hh(6, 6), f(6), df(6, 6)
        hh = 0; df = 0
        call fo_scalar_first( &
            p(1), &
            p(2), &
            p(3), &
            mass, charge, e0, &
            aj(1, 1), &
            aj(1, 2), &
            aj(1, 3), &
            aj(1, 4), &
            aj(1, 5), &
            aj(1, 6), &
            aj(1, 7), &
            aj(1, 8), &
            aj(1, 9), &
            aj(1, 10), &
            aj(2, 1), &
            aj(2, 2), &
            aj(2, 3), &
            aj(2, 4), &
            aj(2, 5), &
            aj(2, 6), &
            aj(2, 7), &
            aj(2, 8), &
            aj(2, 9), &
            aj(2, 10), &
            aj(3, 1), &
            aj(3, 2), &
            aj(3, 3), &
            aj(3, 4), &
            aj(3, 5), &
            aj(3, 6), &
            aj(3, 7), &
            aj(3, 8), &
            aj(3, 9), &
            aj(3, 10), &
            gj(1, 1), &
            gj(1, 2), &
            gj(1, 3), &
            gj(1, 4), &
            gj(1, 5), &
            gj(1, 6), &
            gj(1, 7), &
            gj(1, 8), &
            gj(1, 9), &
            gj(1, 10), &
            gj(2, 1), &
            gj(2, 2), &
            gj(2, 3), &
            gj(2, 4), &
            gj(2, 5), &
            gj(2, 6), &
            gj(2, 7), &
            gj(2, 8), &
            gj(2, 9), &
            gj(2, 10), &
            gj(3, 1), &
            gj(3, 2), &
            gj(3, 3), &
            gj(3, 4), &
            gj(3, 5), &
            gj(3, 6), &
            gj(3, 7), &
            gj(3, 8), &
            gj(3, 9), &
            gj(3, 10), &
            gj(4, 1), &
            gj(4, 2), &
            gj(4, 3), &
            gj(4, 4), &
            gj(4, 5), &
            gj(4, 6), &
            gj(4, 7), &
            gj(4, 8), &
            gj(4, 9), &
            gj(4, 10), &
            gj(5, 1), &
            gj(5, 2), &
            gj(5, 3), &
            gj(5, 4), &
            gj(5, 5), &
            gj(5, 6), &
            gj(5, 7), &
            gj(5, 8), &
            gj(5, 9), &
            gj(5, 10), &
            gj(6, 1), &
            gj(6, 2), &
            gj(6, 3), &
            gj(6, 4), &
            gj(6, 5), &
            gj(6, 6), &
            gj(6, 7), &
            gj(6, 8), &
            gj(6, 9), &
            gj(6, 10), &
            tj(1), &
            tj(2), &
            tj(3), &
            tj(4), &
            tj(5), &
            tj(6), &
            tj(7), &
            tj(8), &
            tj(9), &
            tj(10), &
            energy, kenergy, &
            gh(1), &
            gh(2), &
            gh(3), &
            gh(4), &
            gh(5), &
            gh(6), &
            f(1), &
            f(2), &
            f(3), &
            f(4), &
            f(5), &
            f(6))
    end subroutine fo_eval_first
subroutine fo_radial(gj,rj)
real(dp),intent(in)::gj(6,10)
real(dp),intent(out)::rj(10)
call fo_radial_scalar( &
gj(1,1), &
gj(1,2), &
gj(1,3), &
gj(1,4), &
gj(1,5), &
gj(1,6), &
gj(1,7), &
gj(1,8), &
gj(1,9), &
gj(1,10), &
gj(2,1), &
gj(2,2), &
gj(2,3), &
gj(2,4), &
gj(2,5), &
gj(2,6), &
gj(2,7), &
gj(2,8), &
gj(2,9), &
gj(2,10), &
gj(3,1), &
gj(3,2), &
gj(3,3), &
gj(3,4), &
gj(3,5), &
gj(3,6), &
gj(3,7), &
gj(3,8), &
gj(3,9), &
gj(3,10), &
gj(4,1), &
gj(4,2), &
gj(4,3), &
gj(4,4), &
gj(4,5), &
gj(4,6), &
gj(4,7), &
gj(4,8), &
gj(4,9), &
gj(4,10), &
gj(5,1), &
gj(5,2), &
gj(5,3), &
gj(5,4), &
gj(5,5), &
gj(5,6), &
gj(5,7), &
gj(5,8), &
gj(5,9), &
gj(5,10), &
gj(6,1), &
gj(6,2), &
gj(6,3), &
gj(6,4), &
gj(6,5), &
gj(6,6), &
gj(6,7), &
gj(6,8), &
gj(6,9), &
gj(6,10), &
rj(1), &
rj(2), &
rj(3), &
rj(4), &
rj(5), &
rj(6), &
rj(7), &
rj(8), &
rj(9), &
rj(10))
end subroutine fo_radial
subroutine fo_metric_monitor(rj,tj,mref,out)
real(dp),intent(in)::rj(10),tj(3),mref
real(dp),intent(out)::out(10)
call fo_metric_monitor_scalar( &
rj(1), &
rj(2), &
rj(3), &
rj(4), &
rj(5), &
rj(6), &
rj(7), &
rj(8), &
rj(9), &
rj(10), &
tj(1),tj(2),tj(3),mref, &
out(1), &
out(2), &
out(3), &
out(4), &
out(5), &
out(6), &
out(7), &
out(8), &
out(9), &
out(10))
end subroutine fo_metric_monitor
include 'fo_scalar.inc'
    include 'fo_scalar_first.inc'
include 'fo_monitor.inc'
include 'fo_product.inc'
include 'fo_radial.inc'
include 'fo_metric_monitor.inc'
end module fo_kernel
