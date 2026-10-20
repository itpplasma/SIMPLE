module fo_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    subroutine &
        fo_eval(p, aj, gj, mass, charge, e0, tj, energy, kenergy, gh, hh, f, df)
        real(dp), intent(in)::p(3), aj(3, 10), gj(6, 10), mass, charge, e0, tj(10)
        real(dp), intent(out)::energy, kenergy, gh(6), hh(6, 6), f(6), df(6, 6)
        call fo_scalar(p(1), p(2), p(3), mass, charge, e0, aj(1, 1), &
                 aj(1, 2), aj(1, 3), aj(1, 4), aj(1, 5), aj(1, 6), aj(1, 7), aj(1, 8), &
                aj(1, 9), aj(1, 10), aj(2, 1), aj(2, 2), aj(2, 3), aj(2, 4), aj(2, 5), &
                aj(2, 6), aj(2, 7), aj(2, 8), aj(2, 9), aj(2, 10), aj(3, 1), aj(3, 2), &
                 aj(3, 3), aj(3, 4), aj(3, 5), aj(3, 6), aj(3, 7), aj(3, 8), aj(3, 9), &
                aj(3, 10), gj(1, 1), gj(1, 2), gj(1, 3), gj(1, 4), gj(1, 5), gj(1, 6), &
                gj(1, 7), gj(1, 8), gj(1, 9), gj(1, 10), gj(2, 1), gj(2, 2), gj(2, 3), &
                gj(2, 4), gj(2, 5), gj(2, 6), gj(2, 7), gj(2, 8), gj(2, 9), gj(2, 10), &
                 gj(3, 1), gj(3, 2), gj(3, 3), gj(3, 4), gj(3, 5), gj(3, 6), gj(3, 7), &
                gj(3, 8), gj(3, 9), gj(3, 10), gj(4, 1), gj(4, 2), gj(4, 3), gj(4, 4), &
                gj(4, 5), gj(4, 6), gj(4, 7), gj(4, 8), gj(4, 9), gj(4, 10), gj(5, 1), &
                 gj(5, 2), gj(5, 3), gj(5, 4), gj(5, 5), gj(5, 6), gj(5, 7), gj(5, 8), &
                gj(5, 9), gj(5, 10), gj(6, 1), gj(6, 2), gj(6, 3), gj(6, 4), gj(6, 5), &
               gj(6, 6), gj(6, 7), gj(6, 8), gj(6, 9), gj(6, 10), tj(1), tj(2), tj(3), &
                    tj(4), tj(5), tj(6), tj(7), tj(8), tj(9), tj(10), energy, kenergy, &
               gh(1), gh(2), gh(3), gh(4), gh(5), gh(6), hh(1, 1), hh(2, 1), hh(3, 1), &
                 hh(4, 1), hh(5, 1), hh(6, 1), hh(1, 2), hh(2, 2), hh(3, 2), hh(4, 2), &
                 hh(5, 2), hh(6, 2), hh(1, 3), hh(2, 3), hh(3, 3), hh(4, 3), hh(5, 3), &
                 hh(6, 3), hh(1, 4), hh(2, 4), hh(3, 4), hh(4, 4), hh(5, 4), hh(6, 4), &
                 hh(1, 5), hh(2, 5), hh(3, 5), hh(4, 5), hh(5, 5), hh(6, 5), hh(1, 6), &
             hh(2, 6), hh(3, 6), hh(4, 6), hh(5, 6), hh(6, 6), f(1), f(2), f(3), f(4), &
               f(5), f(6), df(1, 1), df(2, 1), df(3, 1), df(4, 1), df(5, 1), df(6, 1), &
                 df(1, 2), df(2, 2), df(3, 2), df(4, 2), df(5, 2), df(6, 2), df(1, 3), &
                 df(2, 3), df(3, 3), df(4, 3), df(5, 3), df(6, 3), df(1, 4), df(2, 4), &
                 df(3, 4), df(4, 4), df(5, 4), df(6, 4), df(1, 5), df(2, 5), df(3, 5), &
                 df(4, 5), df(5, 5), df(6, 5), df(1, 6), df(2, 6), df(3, 6), df(4, 6), &
                       df(5, 6), df(6, 6))
    end subroutine fo_eval
    subroutine &
        fo_eval_first(p, aj, gj, mass, charge, e0, tj, energy, kenergy, gh, hh, f, df)
        real(dp), intent(in)::p(3), aj(3, 10), gj(6, 10), mass, charge, e0, tj(10)
        real(dp), intent(out)::energy, kenergy, gh(6), hh(6, 6), f(6), df(6, 6)
        real(dp) :: work1, work2, work3, work4, work5, work6, work7, work8, &
                    work9, work10, work11, work12, work13, work14, work15, work16, &
                    work17, work18, work19, work20, work21, work22, work23, work24, &
                    work25, work26, work27, work28, work29

        hh = 0; df = 0

        work1 = gj(5, 1)**2
        work2 = gj(3, 1)**2
        work3 = gj(2, 1)**2
        work4 = gj(1, 1)*gj(4, 1)*gj(6, 1) - gj(1, 1)*work1 + &
                gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*work2 - gj(6, 1)*work3
        work5 = p(1) - aj(1, 1)*charge
        work6 = work5**2
        work7 = gj(4, 1)*gj(6, 1) - work1
        work8 = work6*work7*1.0_dp/2.0_dp/(mass*work4)
        work9 = p(2) - aj(2, 1)*charge
        work10 = work9**2
        work11 = gj(1, 1)*gj(6, 1) - work2
        work12 = work10*work11*1.0_dp/2.0_dp/(mass*work4)
        work13 = p(3) - aj(3, 1)*charge
        work14 = work13**2
        work15 = gj(1, 1)*gj(4, 1) - work3
        work16 = work14*work15*1.0_dp/2.0_dp/(mass*work4)
        work17 = -gj(2, 1)*gj(6, 1) + gj(3, 1)*gj(5, 1)
        work18 = work5*work9*work17/(mass*work4)
        work19 = gj(2, 1)*gj(5, 1) - gj(3, 1)*gj(4, 1)
        work20 = work5*work13*work19/(mass*work4)
        work21 = -gj(1, 1)*gj(5, 1) + gj(2, 1)*gj(3, 1)
        work22 = work9*work13*work21/(mass*work4)
        work23 = -e0 + work8 + work12 + work16 + work18 + work20 + work22
        work24 = gj(1, 1)*gj(4, 1)*gj(6, 2) + gj(1, 1)*gj(4, 2)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)*gj(5, 2)*2 + gj(1, 2)*gj(4, 1)*gj(6, 1) - gj(1, 2)*work1 &
                 - gj(2, 1)*gj(2, 2)*gj(6, 1)*2 + gj(2, 1)*gj(3, 1)*gj(5, 2)*2 + &
                 gj(2, 1)*gj(3, 2)*gj(5, 1)*2 + gj(2, 2)*gj(3, 1)*gj(5, 1)*2 - &
                 gj(3, 1)*gj(3, 2)*gj(4, 1)*2 - gj(4, 2)*work2 - gj(6, 2)*work3
        work25 = -aj(1, 2)*charge*work5*work7/(mass*work4) - &
                 aj(1, 2)*charge*work9*work17/(mass*work4) - aj(1, 2)* &
                 charge*work13*work19/(mass*work4) - &
                 aj(2, 2)*charge*work5*work17/(mass*work4) - aj(2, 2)*charge*work9* &
                 work11/(mass*work4) - aj(2, 2)*charge*work13*work21/(mass*work4) - &
                 aj(3, 2)*charge*work5*work19/(mass*work4) - &
                 aj(3, 2)*charge*work9*work21/(mass*work4) - &
                 aj(3, 2)*charge*work13*work15/(mass*work4) - work6*work24* &
                 work7*1.0_dp/2.0_dp/(mass*work4**2) - &
                 work10*work24*work11*1.0_dp/2.0_dp/(mass*work4**2) - work14* &
                 work15*work24*1.0_dp/2.0_dp/(mass*work4**2) - &
                 work5*work9*work24*work17/(mass*work4**2) - work5* &
                 work13*work24*work19/(mass*work4**2) - &
                 work9*work13*work24*work21/(mass*work4**2) + work6*(gj(4, 1)*gj(6, 2) &
               + gj(4, 2)*gj(6, 1) - gj(5, 1)*gj(5, 2)*2)*1.0_dp/2.0_dp/(mass*work4) + &
                 work10*(gj(1, 1)*gj(6, 2) + gj(1, 2)*gj(6, 1) - gj(3, 1)* &
                  gj(3, 2)*2)*1.0_dp/2.0_dp/(mass*work4) + work14*(gj(1, 1)*gj(4, 2) + &
                 gj(1, 2)*gj(4, 1) - gj(2, 1)*gj(2, 2)*2)*1.0_dp/2.0_dp/(mass*work4) + &
             work5*work9*(-gj(2, 1)*gj(6, 2) - gj(2, 2)*gj(6, 1) + gj(3, 1)*gj(5, 2) + &
                   gj(3, 2)*gj(5, 1))/(mass*work4) + work5*work13*(gj(2, 1)*gj(5, 2) + &
             gj(2, 2)*gj(5, 1) - gj(3, 1)*gj(4, 2) - gj(3, 2)*gj(4, 1))/(mass*work4) + &
            work9*work13*(-gj(1, 1)*gj(5, 2) - gj(1, 2)*gj(5, 1) + gj(2, 1)*gj(3, 2) + &
                               gj(2, 2)*gj(3, 1))/(mass*work4)
        work26 = gj(1, 1)*gj(4, 1)*gj(6, 3) + gj(1, 1)*gj(4, 3)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)*gj(5, 3)*2 + gj(1, 3)*gj(4, 1)*gj(6, 1) - gj(1, 3)*work1 &
                 - gj(2, 1)*gj(2, 3)*gj(6, 1)*2 + gj(2, 1)*gj(3, 1)*gj(5, 3)*2 + &
                 gj(2, 1)*gj(3, 3)*gj(5, 1)*2 + gj(2, 3)*gj(3, 1)*gj(5, 1)*2 - &
                 gj(3, 1)*gj(3, 3)*gj(4, 1)*2 - gj(4, 3)*work2 - gj(6, 3)*work3
        work27 = -aj(1, 3)*charge*work5*work7/(mass*work4) - &
                 aj(1, 3)*charge*work9*work17/(mass*work4) - aj(1, 3)* &
                 charge*work13*work19/(mass*work4) - &
                 aj(2, 3)*charge*work5*work17/(mass*work4) - aj(2, 3)*charge*work9* &
                 work11/(mass*work4) - aj(2, 3)*charge*work13*work21/(mass*work4) - &
                 aj(3, 3)*charge*work5*work19/(mass*work4) - &
                 aj(3, 3)*charge*work9*work21/(mass*work4) - &
                 aj(3, 3)*charge*work13*work15/(mass*work4) - work6*work26* &
                 work7*1.0_dp/2.0_dp/(mass*work4**2) - &
                 work10*work26*work11*1.0_dp/2.0_dp/(mass*work4**2) - work14* &
                 work15*work26*1.0_dp/2.0_dp/(mass*work4**2) - &
                 work5*work9*work26*work17/(mass*work4**2) - work5* &
                 work13*work26*work19/(mass*work4**2) - &
                 work9*work13*work26*work21/(mass*work4**2) + work6*(gj(4, 1)*gj(6, 3) &
               + gj(4, 3)*gj(6, 1) - gj(5, 1)*gj(5, 3)*2)*1.0_dp/2.0_dp/(mass*work4) + &
                 work10*(gj(1, 1)*gj(6, 3) + gj(1, 3)*gj(6, 1) - gj(3, 1)* &
                  gj(3, 3)*2)*1.0_dp/2.0_dp/(mass*work4) + work14*(gj(1, 1)*gj(4, 3) + &
                 gj(1, 3)*gj(4, 1) - gj(2, 1)*gj(2, 3)*2)*1.0_dp/2.0_dp/(mass*work4) + &
             work5*work9*(-gj(2, 1)*gj(6, 3) - gj(2, 3)*gj(6, 1) + gj(3, 1)*gj(5, 3) + &
                   gj(3, 3)*gj(5, 1))/(mass*work4) + work5*work13*(gj(2, 1)*gj(5, 3) + &
             gj(2, 3)*gj(5, 1) - gj(3, 1)*gj(4, 3) - gj(3, 3)*gj(4, 1))/(mass*work4) + &
            work9*work13*(-gj(1, 1)*gj(5, 3) - gj(1, 3)*gj(5, 1) + gj(2, 1)*gj(3, 3) + &
                               gj(2, 3)*gj(3, 1))/(mass*work4)
        work28 = gj(1, 1)*gj(4, 1)*gj(6, 4) + gj(1, 1)*gj(4, 4)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)*gj(5, 4)*2 + gj(1, 4)*gj(4, 1)*gj(6, 1) - gj(1, 4)*work1 &
                 - gj(2, 1)*gj(2, 4)*gj(6, 1)*2 + gj(2, 1)*gj(3, 1)*gj(5, 4)*2 + &
                 gj(2, 1)*gj(3, 4)*gj(5, 1)*2 + gj(2, 4)*gj(3, 1)*gj(5, 1)*2 - &
                 gj(3, 1)*gj(3, 4)*gj(4, 1)*2 - gj(4, 4)*work2 - gj(6, 4)*work3
        work29 = -aj(1, 4)*charge*work5*work7/(mass*work4) - &
                 aj(1, 4)*charge*work9*work17/(mass*work4) - aj(1, 4)* &
                 charge*work13*work19/(mass*work4) - &
                 aj(2, 4)*charge*work5*work17/(mass*work4) - aj(2, 4)*charge*work9* &
                 work11/(mass*work4) - aj(2, 4)*charge*work13*work21/(mass*work4) - &
                 aj(3, 4)*charge*work5*work19/(mass*work4) - &
                 aj(3, 4)*charge*work9*work21/(mass*work4) - &
                 aj(3, 4)*charge*work13*work15/(mass*work4) - work6*work28* &
                 work7*1.0_dp/2.0_dp/(mass*work4**2) - &
                 work10*work28*work11*1.0_dp/2.0_dp/(mass*work4**2) - work14* &
                 work15*work28*1.0_dp/2.0_dp/(mass*work4**2) - &
                 work5*work9*work28*work17/(mass*work4**2) - work5* &
                 work13*work28*work19/(mass*work4**2) - &
                 work9*work13*work28*work21/(mass*work4**2) + work6*(gj(4, 1)*gj(6, 4) &
               + gj(4, 4)*gj(6, 1) - gj(5, 1)*gj(5, 4)*2)*1.0_dp/2.0_dp/(mass*work4) + &
                 work10*(gj(1, 1)*gj(6, 4) + gj(1, 4)*gj(6, 1) - gj(3, 1)* &
                  gj(3, 4)*2)*1.0_dp/2.0_dp/(mass*work4) + work14*(gj(1, 1)*gj(4, 4) + &
                 gj(1, 4)*gj(4, 1) - gj(2, 1)*gj(2, 4)*2)*1.0_dp/2.0_dp/(mass*work4) + &
             work5*work9*(-gj(2, 1)*gj(6, 4) - gj(2, 4)*gj(6, 1) + gj(3, 1)*gj(5, 4) + &
                   gj(3, 4)*gj(5, 1))/(mass*work4) + work5*work13*(gj(2, 1)*gj(5, 4) + &
             gj(2, 4)*gj(5, 1) - gj(3, 1)*gj(4, 4) - gj(3, 4)*gj(4, 1))/(mass*work4) + &
            work9*work13*(-gj(1, 1)*gj(5, 4) - gj(1, 4)*gj(5, 1) + gj(2, 1)*gj(3, 4) + &
                               gj(2, 4)*gj(3, 1))/(mass*work4)
        energy = work8 + work12 + work16 + work18 + work20 + work22
        kenergy = tj(1)*work23
        gh(1) = work25
        gh(2) = work27
        gh(3) = work29
        gh(4) = (-aj(1, 1)*charge*gj(4, 1)*gj(6, 1) + aj(1, 1)*charge*work1 + &
               aj(2, 1)*charge*gj(2, 1)*gj(6, 1) - aj(2, 1)*charge*gj(3, 1)*gj(5, 1) - &
               aj(3, 1)*charge*gj(2, 1)*gj(5, 1) + aj(3, 1)*charge*gj(3, 1)*gj(4, 1) + &
            gj(2, 1)*gj(5, 1)*p(3) - gj(2, 1)*gj(6, 1)*p(2) - gj(3, 1)*gj(4, 1)*p(3) + &
                 gj(3, 1)*gj(5, 1)*p(2) + gj(4, 1)*gj(6, 1)*p(1) - &
                 p(1)*work1)/(gj(1, 1)*gj(4, 1)*gj(6, 1)*mass - gj(1, 1)*mass*work1 + &
                             gj(2, 1)*gj(3, 1)*gj(5, 1)*mass*2 - gj(4, 1)*mass*work2 - &
                              gj(6, 1)*mass*work3)
        gh(5) = (aj(1, 1)*charge*gj(2, 1)*gj(6, 1) - &
               aj(1, 1)*charge*gj(3, 1)*gj(5, 1) - aj(2, 1)*charge*gj(1, 1)*gj(6, 1) + &
                 aj(2, 1)*charge*work2 + aj(3, 1)*charge*gj(1, 1)*gj(5, 1) - &
                 aj(3, 1)*charge*gj(2, 1)*gj(3, 1) - gj(1, 1)*gj(5, 1)*p(3) + &
            gj(1, 1)*gj(6, 1)*p(2) + gj(2, 1)*gj(3, 1)*p(3) - gj(2, 1)*gj(6, 1)*p(1) + &
               gj(3, 1)*gj(5, 1)*p(1) - p(2)*work2)/(gj(1, 1)*gj(4, 1)*gj(6, 1)*mass - &
                             gj(1, 1)*mass*work1 + gj(2, 1)*gj(3, 1)*gj(5, 1)*mass*2 - &
                                              gj(4, 1)*mass*work2 - gj(6, 1)*mass*work3)
        gh(6) = (-aj(1, 1)*charge*gj(2, 1)*gj(5, 1) + &
               aj(1, 1)*charge*gj(3, 1)*gj(4, 1) + aj(2, 1)*charge*gj(1, 1)*gj(5, 1) - &
               aj(2, 1)*charge*gj(2, 1)*gj(3, 1) - aj(3, 1)*charge*gj(1, 1)*gj(4, 1) + &
             aj(3, 1)*charge*work3 + gj(1, 1)*gj(4, 1)*p(3) - gj(1, 1)*gj(5, 1)*p(2) + &
            gj(2, 1)*gj(3, 1)*p(2) + gj(2, 1)*gj(5, 1)*p(1) - gj(3, 1)*gj(4, 1)*p(1) - &
                 p(3)*work3)/(gj(1, 1)*gj(4, 1)*gj(6, 1)*mass - gj(1, 1)*mass*work1 + &
                             gj(2, 1)*gj(3, 1)*gj(5, 1)*mass*2 - gj(4, 1)*mass*work2 - &
                              gj(6, 1)*mass*work3)
        f(1) = (-aj(1, 1)*charge*gj(4, 1)*gj(6, 1)*tj(1) + &
               aj(1, 1)*charge*tj(1)*work1 + aj(2, 1)*charge*gj(2, 1)*gj(6, 1)*tj(1) - &
                aj(2, 1)*charge*gj(3, 1)*gj(5, 1)*tj(1) - &
                aj(3, 1)*charge*gj(2, 1)*gj(5, 1)*tj(1) + &
              aj(3, 1)*charge*gj(3, 1)*gj(4, 1)*tj(1) + gj(2, 1)*gj(5, 1)*p(3)*tj(1) - &
                gj(2, 1)*gj(6, 1)*p(2)*tj(1) - gj(3, 1)*gj(4, 1)*p(3)*tj(1) + &
                gj(3, 1)*gj(5, 1)*p(2)*tj(1) + gj(4, 1)*gj(6, 1)*p(1)*tj(1) - &
                p(1)*tj(1)*work1)/(gj(1, 1)*gj(4, 1)*gj(6, 1)*mass - &
                             gj(1, 1)*mass*work1 + gj(2, 1)*gj(3, 1)*gj(5, 1)*mass*2 - &
                                   gj(4, 1)*mass*work2 - gj(6, 1)*mass*work3)
        f(2) = (aj(1, 1)*charge*gj(2, 1)*gj(6, 1)*tj(1) - &
                aj(1, 1)*charge*gj(3, 1)*gj(5, 1)*tj(1) - &
               aj(2, 1)*charge*gj(1, 1)*gj(6, 1)*tj(1) + aj(2, 1)*charge*tj(1)*work2 + &
                aj(3, 1)*charge*gj(1, 1)*gj(5, 1)*tj(1) - &
              aj(3, 1)*charge*gj(2, 1)*gj(3, 1)*tj(1) - gj(1, 1)*gj(5, 1)*p(3)*tj(1) + &
                gj(1, 1)*gj(6, 1)*p(2)*tj(1) + gj(2, 1)*gj(3, 1)*p(3)*tj(1) - &
                gj(2, 1)*gj(6, 1)*p(1)*tj(1) + gj(3, 1)*gj(5, 1)*p(1)*tj(1) - &
                p(2)*tj(1)*work2)/(gj(1, 1)*gj(4, 1)*gj(6, 1)*mass - &
                             gj(1, 1)*mass*work1 + gj(2, 1)*gj(3, 1)*gj(5, 1)*mass*2 - &
                                   gj(4, 1)*mass*work2 - gj(6, 1)*mass*work3)
        f(3) = (-aj(1, 1)*charge*gj(2, 1)*gj(5, 1)*tj(1) + &
                aj(1, 1)*charge*gj(3, 1)*gj(4, 1)*tj(1) + &
                aj(2, 1)*charge*gj(1, 1)*gj(5, 1)*tj(1) - &
                aj(2, 1)*charge*gj(2, 1)*gj(3, 1)*tj(1) - &
               aj(3, 1)*charge*gj(1, 1)*gj(4, 1)*tj(1) + aj(3, 1)*charge*tj(1)*work3 + &
                gj(1, 1)*gj(4, 1)*p(3)*tj(1) - gj(1, 1)*gj(5, 1)*p(2)*tj(1) + &
                gj(2, 1)*gj(3, 1)*p(2)*tj(1) + gj(2, 1)*gj(5, 1)*p(1)*tj(1) - &
                gj(3, 1)*gj(4, 1)*p(1)*tj(1) - p(3)*tj(1)* &
                work3)/(gj(1, 1)*gj(4, 1)*gj(6, 1)*mass - gj(1, 1)*mass*work1 + &
              gj(2, 1)*gj(3, 1)*gj(5, 1)*mass*2 - gj(4, 1)*mass*work2 - gj(6, 1)*mass* &
                        work3)
        f(4) = -(tj(1)*work25 + tj(2)*work23)
        f(5) = -(tj(1)*work27 + tj(3)*work23)
        f(6) = -(tj(1)*work29 + tj(4)*work23)

    end subroutine fo_eval_first
    subroutine fo_radial(gj, rj)
        real(dp), intent(in)::gj(6, 10)
        real(dp), intent(out)::rj(10)
        call fo_radial_scalar(gj(1, 1), gj(1, 2), gj(1, 3), gj(1, 4), gj(1, 5), &
                gj(1, 6), gj(1, 7), gj(1, 8), gj(1, 9), gj(1, 10), gj(2, 1), gj(2, 2), &
                 gj(2, 3), gj(2, 4), gj(2, 5), gj(2, 6), gj(2, 7), gj(2, 8), gj(2, 9), &
                gj(2, 10), gj(3, 1), gj(3, 2), gj(3, 3), gj(3, 4), gj(3, 5), gj(3, 6), &
                gj(3, 7), gj(3, 8), gj(3, 9), gj(3, 10), gj(4, 1), gj(4, 2), gj(4, 3), &
                gj(4, 4), gj(4, 5), gj(4, 6), gj(4, 7), gj(4, 8), gj(4, 9), gj(4, 10), &
                 gj(5, 1), gj(5, 2), gj(5, 3), gj(5, 4), gj(5, 5), gj(5, 6), gj(5, 7), &
                gj(5, 8), gj(5, 9), gj(5, 10), gj(6, 1), gj(6, 2), gj(6, 3), gj(6, 4), &
            gj(6, 5), gj(6, 6), gj(6, 7), gj(6, 8), gj(6, 9), gj(6, 10), rj(1), rj(2), &
                              rj(3), rj(4), rj(5), rj(6), rj(7), rj(8), rj(9), rj(10))
    end subroutine fo_radial
    subroutine fo_metric_monitor(rj, tj, mref, out)
        real(dp), intent(in)::rj(10), tj(3), mref
        real(dp), intent(out)::out(10)
        call fo_metric_monitor_scalar(rj(1), rj(2), rj(3), rj(4), rj(5), &
                rj(6), rj(7), rj(8), rj(9), rj(10), tj(1), tj(2), tj(3), mref, out(1), &
                       out(2), out(3), out(4), out(5), out(6), out(7), out(8), out(9), &
                                      out(10))
    end subroutine fo_metric_monitor
    include 'fo_scalar.inc'
    include 'fo_scalar_first.inc'
    include 'fo_monitor.inc'
    include 'fo_product.inc'
    include 'fo_radial.inc'
    include 'fo_metric_monitor.inc'
    subroutine &
       fo_eval_factored(p, aj, gj, mass, charge, e0, tj, energy, kenergy, gh, hh, f, df)
        real(dp), intent(in)::p(3), aj(3, 10), gj(6, 10), mass, charge, e0, tj(10)
        real(dp), intent(out)::energy, kenergy, gh(6), hh(6, 6), f(6), df(6, 6)
        real(dp)::inverse(6), velocity(3)
        call &
          fo_eval_first(p, aj, gj, mass, charge, e0, tj, energy, kenergy, gh, hh, f, df)
        block
            real(dp) :: work1, work2, work3, work4, work5, work6, work7, work8, &
                        work9

            work1 = gj(4, 1)*gj(6, 1) - gj(5, 1)*gj(5, 1)
            work2 = -gj(2, 1)*gj(6, 1) + gj(3, 1)*gj(5, 1)
            work3 = gj(2, 1)*gj(5, 1) - gj(3, 1)*gj(4, 1)
            work4 = gj(1, 1)*gj(6, 1) - gj(3, 1)*gj(3, 1)
            work5 = -gj(1, 1)*gj(5, 1) + gj(2, 1)*gj(3, 1)
            work6 = gj(1, 1)*gj(4, 1) - gj(2, 1)*gj(2, 1)
            work7 = p(1) - aj(1, 1)*charge
            work8 = p(2) - aj(2, 1)*charge
            work9 = p(3) - aj(3, 1)*charge
            inverse(1) = work1/(gj(1, 1)*gj(4, 1)*gj(6, 1) - gj(1, 1)*gj(5, 1)**2 + &
                       gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 - gj(6, 1)* &
                                gj(2, 1)**2)
            inverse(2) = work2/(gj(1, 1)*gj(4, 1)*gj(6, 1) - gj(1, 1)*gj(5, 1)**2 + &
                       gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 - gj(6, 1)* &
                                gj(2, 1)**2)
            inverse(3) = work3/(gj(1, 1)*gj(4, 1)*gj(6, 1) - gj(1, 1)*gj(5, 1)**2 + &
                       gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 - gj(6, 1)* &
                                gj(2, 1)**2)
            inverse(4) = work4/(gj(1, 1)*gj(4, 1)*gj(6, 1) - gj(1, 1)*gj(5, 1)**2 + &
                       gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 - gj(6, 1)* &
                                gj(2, 1)**2)
            inverse(5) = work5/(gj(1, 1)*gj(4, 1)*gj(6, 1) - gj(1, 1)*gj(5, 1)**2 + &
                       gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 - gj(6, 1)* &
                                gj(2, 1)**2)
            inverse(6) = work6/(gj(1, 1)*gj(4, 1)*gj(6, 1) - gj(1, 1)*gj(5, 1)**2 + &
                       gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 - gj(6, 1)* &
                                gj(2, 1)**2)
            velocity(1) = work7*work1/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                   - gj(6, 1)*gj(2, 1)**2) + work8*work2/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                   - gj(6, 1)*gj(2, 1)**2) + work9*work3/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                                                             - gj(6, 1)*gj(2, 1)**2) + 0
            velocity(2) = work7*work2/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                   - gj(6, 1)*gj(2, 1)**2) + work8*work4/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                   - gj(6, 1)*gj(2, 1)**2) + work9*work5/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                                                             - gj(6, 1)*gj(2, 1)**2) + 0
            velocity(3) = work7*work3/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                   - gj(6, 1)*gj(2, 1)**2) + work8*work5/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                   - gj(6, 1)*gj(2, 1)**2) + work9*work6/(gj(1, 1)*gj(4, 1)*gj(6, 1) - &
            gj(1, 1)*gj(5, 1)**2 + gj(2, 1)*gj(3, 1)*gj(5, 1)*2 - gj(4, 1)*gj(3, 1)**2 &
                                                             - gj(6, 1)*gj(2, 1)**2) + 0

        end block
        block
            real(dp) :: work1, work2, work3, work4, work5, work6, work7, work8, &
                        work9, work10, work11, work12, work13, work14, work15, work16, &
                       work17, work18, work19, work20, work21, work22, work23, work24, &
                       work25, work26, work27, work28, work29, work30, work31, work32, &
                        work33

            work1 = aj(1, 2)*charge + gj(1, 2)*velocity(1) + gj(2, 2)*velocity(2) + &
                    gj(3, 2)*velocity(3)
            work2 = aj(2, 2)*charge + gj(2, 2)*velocity(1) + gj(4, 2)*velocity(2) + &
                    gj(5, 2)*velocity(3)
            work3 = inverse(2)*work1*work2
            work4 = aj(3, 2)*charge + gj(3, 2)*velocity(1) + gj(5, 2)*velocity(2) + &
                    gj(6, 2)*velocity(3)
            work5 = inverse(3)*work1*work4
            work6 = inverse(5)*work2*work4
            work7 = aj(1, 3)*charge + gj(1, 3)*velocity(1) + gj(2, 3)*velocity(2) + &
                    gj(3, 3)*velocity(3)
            work8 = aj(2, 3)*charge + gj(2, 3)*velocity(1) + gj(4, 3)*velocity(2) + &
                    gj(5, 3)*velocity(3)
            work9 = aj(3, 3)*charge + gj(3, 3)*velocity(1) + gj(5, 3)*velocity(2) + &
                    gj(6, 3)*velocity(3)
            work10 = (-aj(1, 6)*charge*velocity(1) - aj(2, 6)*charge*velocity(2) - &
                      aj(3, 6)*charge*velocity(3) - gj(1, 6)* &
                      velocity(1)*velocity(1)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(2, 6)*velocity(1)*velocity(2)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(3, 6)*velocity(1)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(4, 6)*velocity(2)*velocity(2)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(5, 6)*velocity(2)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(6, 6)*velocity(3)*velocity(3)*1/real(2, kind=kind(0.0_dp)) + &
                      inverse(1)*work1*work7 + inverse(2)*work1*work8 + &
                      inverse(2)*work7*work2 + inverse(3)*work1*work9 + &
                      inverse(3)*work7*work4 + inverse(4)*work2*work8 + &
                      inverse(5)*work2*work9 + inverse(5)*work8*work4 + &
                      inverse(6)*work4*work9 + 0)/mass
            work11 = aj(1, 4)*charge + gj(1, 4)*velocity(1) + gj(2, 4)*velocity(2) &
                     + gj(3, 4)*velocity(3)
            work12 = aj(2, 4)*charge + gj(2, 4)*velocity(1) + gj(4, 4)*velocity(2) &
                     + gj(5, 4)*velocity(3)
            work13 = aj(3, 4)*charge + gj(3, 4)*velocity(1) + gj(5, 4)*velocity(2) &
                     + gj(6, 4)*velocity(3)
            work14 = (-aj(1, 8)*charge*velocity(1) - aj(2, 8)*charge*velocity(2) - &
                      aj(3, 8)*charge*velocity(3) - gj(1, 8)* &
                      velocity(1)*velocity(1)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(2, 8)*velocity(1)*velocity(2)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(3, 8)*velocity(1)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(4, 8)*velocity(2)*velocity(2)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(5, 8)*velocity(2)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(6, 8)*velocity(3)*velocity(3)*1/real(2, kind=kind(0.0_dp)) + &
                      inverse(1)*work1*work11 + inverse(2)*work1*work12 + &
                      inverse(2)*work11*work2 + inverse(3)*work1*work13 + &
                      inverse(3)*work11*work4 + inverse(4)*work2*work12 + &
                      inverse(5)*work2*work13 + inverse(5)*work12*work4 + &
                      inverse(6)*work4*work13 + 0)/mass
            work15 = -inverse(1)*work1/mass - inverse(2)*work2/mass - &
                     inverse(3)*work4/mass + 0
            work16 = -inverse(2)*work1/mass - inverse(4)*work2/mass - &
                     inverse(5)*work4/mass + 0
            work17 = -inverse(3)*work1/mass - inverse(5)*work2/mass - &
                     inverse(6)*work4/mass + 0
            work18 = inverse(2)*work7*work8
            work19 = inverse(3)*work7*work9
            work20 = inverse(5)*work8*work9
            work21 = (-aj(1, 9)*charge*velocity(1) - aj(2, 9)*charge*velocity(2) - &
                      aj(3, 9)*charge*velocity(3) - gj(1, 9)* &
                      velocity(1)*velocity(1)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(2, 9)*velocity(1)*velocity(2)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(3, 9)*velocity(1)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(4, 9)*velocity(2)*velocity(2)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(5, 9)*velocity(2)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(6, 9)*velocity(3)*velocity(3)*1/real(2, kind=kind(0.0_dp)) + &
                      inverse(1)*work7*work11 + inverse(2)*work7*work12 + &
                      inverse(2)*work11*work8 + inverse(3)*work7*work13 + &
                      inverse(3)*work11*work9 + inverse(4)*work8*work12 + &
                      inverse(5)*work8*work13 + inverse(5)*work12*work9 + &
                      inverse(6)*work9*work13 + 0)/mass
            work22 = -inverse(1)*work7/mass - inverse(2)*work8/mass - &
                     inverse(3)*work9/mass + 0
            work23 = -inverse(2)*work7/mass - inverse(4)*work8/mass - &
                     inverse(5)*work9/mass + 0
            work24 = -inverse(3)*work7/mass - inverse(5)*work8/mass - &
                     inverse(6)*work9/mass + 0
            work25 = inverse(2)*work11*work12
            work26 = inverse(3)*work11*work13
            work27 = inverse(5)*work12*work13
            work28 = -inverse(1)*work11/mass - inverse(2)*work12/mass - &
                     inverse(3)*work13/mass + 0
            work29 = -inverse(2)*work11/mass - inverse(4)*work12/mass - &
                     inverse(5)*work13/mass + 0
            work30 = -inverse(3)*work11/mass - inverse(5)*work12/mass - &
                     inverse(6)*work13/mass + 0
            work31 = inverse(2)/mass
            work32 = inverse(3)/mass
            work33 = inverse(5)/mass
            hh(1, 1) = (-aj(1, 5)*charge*velocity(1) - aj(2, 5)*charge*velocity(2) &
                        - aj(3, 5)*charge*velocity(3) - gj(1, 5)* &
                        velocity(1)*velocity(1)*1/real(2, kind=kind(0.0_dp)) - &
                       gj(2, 5)*velocity(1)*velocity(2)*2/real(2, kind=kind(0.0_dp)) - &
                       gj(3, 5)*velocity(1)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                       gj(4, 5)*velocity(2)*velocity(2)*1/real(2, kind=kind(0.0_dp)) - &
                       gj(5, 5)*velocity(2)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                       gj(6, 5)*velocity(3)*velocity(3)*1/real(2, kind=kind(0.0_dp)) + &
                        inverse(1)*work1*work1 + work3 + work3 + work5 + work5 + &
                     inverse(4)*work2*work2 + work6 + work6 + inverse(6)*work4*work4 + &
                        0)/mass
            hh(2, 1) = work10
            hh(3, 1) = work14
            hh(4, 1) = work15
            hh(5, 1) = work16
            hh(6, 1) = work17
            hh(1, 2) = work10
            hh(2, 2) = (-aj(1, 7)*charge*velocity(1) - aj(2, 7)*charge*velocity(2) &
                        - aj(3, 7)*charge*velocity(3) - gj(1, 7)* &
                        velocity(1)*velocity(1)*1/real(2, kind=kind(0.0_dp)) - &
                       gj(2, 7)*velocity(1)*velocity(2)*2/real(2, kind=kind(0.0_dp)) - &
                       gj(3, 7)*velocity(1)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                       gj(4, 7)*velocity(2)*velocity(2)*1/real(2, kind=kind(0.0_dp)) - &
                       gj(5, 7)*velocity(2)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                       gj(6, 7)*velocity(3)*velocity(3)*1/real(2, kind=kind(0.0_dp)) + &
                        inverse(1)*work7*work7 + work18 + work18 + work19 + work19 + &
                   inverse(4)*work8*work8 + work20 + work20 + inverse(6)*work9*work9 + &
                        0)/mass
            hh(3, 2) = work21
            hh(4, 2) = work22
            hh(5, 2) = work23
            hh(6, 2) = work24
            hh(1, 3) = work14
            hh(2, 3) = work21
            hh(3, 3) = (-aj(1, 10)*charge*velocity(1) - &
                        aj(2, 10)*charge*velocity(2) - aj(3, 10)*charge*velocity(3) - &
                      gj(1, 10)*velocity(1)*velocity(1)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(2, 10)*velocity(1)*velocity(2)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(3, 10)*velocity(1)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(4, 10)*velocity(2)*velocity(2)*1/real(2, kind=kind(0.0_dp)) - &
                      gj(5, 10)*velocity(2)*velocity(3)*2/real(2, kind=kind(0.0_dp)) - &
                      gj(6, 10)*velocity(3)*velocity(3)*1/real(2, kind=kind(0.0_dp)) + &
                        inverse(1)*work11*work11 + work25 + work25 + work26 + work26 + &
                        inverse(4)*work12*work12 + work27 + work27 + &
                        inverse(6)*work13*work13 + 0)/mass
            hh(4, 3) = work28
            hh(5, 3) = work29
            hh(6, 3) = work30
            hh(1, 4) = work15
            hh(2, 4) = work22
            hh(3, 4) = work28
            hh(4, 4) = inverse(1)/mass
            hh(5, 4) = work31
            hh(6, 4) = work32
            hh(1, 5) = work16
            hh(2, 5) = work23
            hh(3, 5) = work29
            hh(4, 5) = work31
            hh(5, 5) = inverse(4)/mass
            hh(6, 5) = work33
            hh(1, 6) = work17
            hh(2, 6) = work24
            hh(3, 6) = work30
            hh(4, 6) = work32
            hh(5, 6) = work33
            hh(6, 6) = inverse(6)/mass

        end block
        block
            real(dp) :: work1, work2, work3, work4, work5, work6, work7, work8, &
                        work9, work10, work11, work12, work13, work14, work15, work16, &
                        work17, work18, work19, work20, work21, work22

            work1 = gh(4)*tj(2)
            work2 = gh(5)*tj(2)
            work3 = gh(6)*tj(2)
            work4 = gh(1)*tj(2)
            work5 = energy - e0
            work6 = gh(1)*tj(3)
            work7 = gh(2)*tj(2)
            work8 = tj(6)*work5
            work9 = gh(1)*tj(4)
            work10 = gh(3)*tj(2)
            work11 = tj(8)*work5
            work12 = gh(4)*tj(3)
            work13 = gh(5)*tj(3)
            work14 = gh(6)*tj(3)
            work15 = gh(2)*tj(3)
            work16 = gh(2)*tj(4)
            work17 = gh(3)*tj(3)
            work18 = tj(9)*work5
            work19 = gh(4)*tj(4)
            work20 = gh(5)*tj(4)
            work21 = gh(6)*tj(4)
            work22 = gh(3)*tj(4)
            df(1, 1) = work1 + hh(4, 1)*tj(1)
            df(2, 1) = work2 + hh(5, 1)*tj(1)
            df(3, 1) = work3 + hh(6, 1)*tj(1)
            df(4, 1) = -(work4 + work4 + hh(1, 1)*tj(1) + tj(5)*work5)
            df(5, 1) = -(work6 + work7 + hh(2, 1)*tj(1) + work8)
            df(6, 1) = -(work9 + work10 + hh(3, 1)*tj(1) + work11)
            df(1, 2) = work12 + hh(4, 2)*tj(1)
            df(2, 2) = work13 + hh(5, 2)*tj(1)
            df(3, 2) = work14 + hh(6, 2)*tj(1)
            df(4, 2) = -(work6 + work7 + hh(1, 2)*tj(1) + work8)
            df(5, 2) = -(work15 + work15 + hh(2, 2)*tj(1) + tj(7)*work5)
            df(6, 2) = -(work16 + work17 + hh(3, 2)*tj(1) + work18)
            df(1, 3) = work19 + hh(4, 3)*tj(1)
            df(2, 3) = work20 + hh(5, 3)*tj(1)
            df(3, 3) = work21 + hh(6, 3)*tj(1)
            df(4, 3) = -(work9 + work10 + hh(1, 3)*tj(1) + work11)
            df(5, 3) = -(work16 + work17 + hh(2, 3)*tj(1) + work18)
            df(6, 3) = -(work22 + work22 + hh(3, 3)*tj(1) + tj(10)*work5)
            df(1, 4) = hh(4, 4)*tj(1)
            df(2, 4) = hh(5, 4)*tj(1)
            df(3, 4) = hh(6, 4)*tj(1)
            df(4, 4) = -(work1 + hh(1, 4)*tj(1))
            df(5, 4) = -(work12 + hh(2, 4)*tj(1))
            df(6, 4) = -(work19 + hh(3, 4)*tj(1))
            df(1, 5) = hh(4, 5)*tj(1)
            df(2, 5) = hh(5, 5)*tj(1)
            df(3, 5) = hh(6, 5)*tj(1)
            df(4, 5) = -(work2 + hh(1, 5)*tj(1))
            df(5, 5) = -(work13 + hh(2, 5)*tj(1))
            df(6, 5) = -(work20 + hh(3, 5)*tj(1))
            df(1, 6) = hh(4, 6)*tj(1)
            df(2, 6) = hh(5, 6)*tj(1)
            df(3, 6) = hh(6, 6)*tj(1)
            df(4, 6) = -(work3 + hh(1, 6)*tj(1))
            df(5, 6) = -(work14 + hh(2, 6)*tj(1))
            df(6, 6) = -(work21 + hh(3, 6)*tj(1))

        end block
    end subroutine fo_eval_factored
end module fo_kernel
