! Generated array bindings for gen/gen_fo_joint.f90. Do not edit.
module fo_joint_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    subroutine fo_joint_first(p, aj, jets, gj, energy, gh, hh)
        real(dp), intent(in)::p(3), aj(3, 10), jets(0:3, 0:3, 0:3, 2)
        real(dp), intent(out)::gj(6, 10), energy, gh(6), hh(6, 6)
        real(dp) :: work1, work2, work3, work4, work5, work6, work7, &
                    work8, work9, work10, work11, work12, work13, work14, work15, &
                    work16, work17, work18, work19, work20, work21, work22, &
                    work23, work24, work25, work26, work27, work28, work29, &
                    work30, work31, work32, work33, work34, work35, work36, &
                    work37, work38, work39, work40, work41, work42, work43, &
                    work44, work45, work46, work47, work48, work49, work50, &
                    work51, work52, work53, work54, work55, work56, work57, &
                    work58, work59, work60, work61, work62, work63, work64, &
                    work65, work66

        gj = 0; hh = 0

        work1 = p(1) - aj(1, 1)
        work2 = jets(0, 0, 1, 1)*jets(0, 1, 0, 1) + jets(0, 0, 1, 2)* &
                jets(0, 1, 0, 2)
        work3 = jets(1, 0, 0, 1)*jets(1, 0, 0, 1) + jets(1, 0, 0, 2)* &
                jets(1, 0, 0, 2)
        work4 = jets(0, 0, 1, 1)*jets(1, 0, 0, 1) + jets(0, 0, 1, 2)* &
                jets(1, 0, 0, 2)
        work5 = jets(0, 1, 0, 1)*jets(0, 1, 0, 1) + jets(0, 1, 0, 2)* &
                jets(0, 1, 0, 2)
        work6 = jets(0, 1, 0, 1)*jets(1, 0, 0, 1) + jets(0, 1, 0, 2)* &
                jets(1, 0, 0, 2)
        work7 = jets(0, 0, 0, 1)**2 + jets(0, 0, 1, 1)*jets(0, 0, 1, &
                                                  1) + jets(0, 0, 1, 2)*jets(0, 0, 1, 2)
        work8 = p(2) - aj(2, 1)
        work9 = -work7*work6 + work2*work4
        work10 = p(3) - aj(3, 1)
        work11 = work2*work6 - work4*work5
        work12 = work1*(work7*work5 - work2*work2)/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                  work3 + work2*work4*work6*2) + work8*work9/(-work2** &
                                    2*work3 - work4**2*work5 - work6**2*work7 + work7* &
                                   work5*work3 + work2*work4*work6*2) + work10*work11/ &
                 (-work2**2*work3 - work4**2*work5 - work6**2*work7 + &
                  work7*work5*work3 + work2*work4*work6*2) + 0
        work13 = -work2*work3 + work4*work6
        work14 = work1*work9/(-work2**2*work3 - work4**2*work5 &
                              - work6**2*work7 + work7*work5*work3 + work2*work4* &
                              work6*2) + work8*(work7*work3 - work4*work4)/(-work2** &
                                    2*work3 - work4**2*work5 - work6**2*work7 + work7* &
                                   work5*work3 + work2*work4*work6*2) + work10*work13/ &
                 (-work2**2*work3 - work4**2*work5 - work6**2*work7 + &
                  work7*work5*work3 + work2*work4*work6*2) + 0
        work15 = work1*work11/(-work2**2*work3 - work4**2* &
                               work5 - work6**2*work7 + work7*work5*work3 + work2* &
                               work4*work6*2) + work8*work13/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                   work2*work4*work6*2) + work10*(work5*work3 - work6* &
                                  work6)/(-work2**2*work3 - work4**2*work5 - work6**2* &
                                    work7 + work7*work5*work3 + work2*work4*work6*2) + 0
        work16 = jets(0, 0, 1, 1)*jets(1, 0, 1, 1)
        work17 = jets(0, 0, 1, 2)*jets(1, 0, 1, 2)
        work18 = jets(0, 0, 0, 1)*jets(1, 0, 0, 1)*2 + work16 + &
                 work16 + work17 + work17
        work19 = jets(0, 0, 1, 1)*jets(1, 1, 0, 1)
        work20 = jets(0, 0, 1, 2)*jets(1, 1, 0, 2)
        work21 = jets(0, 1, 0, 1)*jets(1, 0, 1, 1)
        work22 = jets(0, 1, 0, 2)*jets(1, 0, 1, 2)
        work23 = work19 + work20 + work21 + work22
        work24 = -work23*work14*work15/real(2, kind=kind(0.0_dp))
        work25 = jets(1, 0, 0, 1)*jets(1, 0, 1, 1)
        work26 = jets(1, 0, 0, 2)*jets(1, 0, 1, 2)
        work27 = jets(0, 0, 1, 1)*jets(2, 0, 0, 1) + jets(0, 0, 1, 2) &
                 *jets(2, 0, 0, 2) + work25 + work26
        work28 = -work27*work12*work15/real(2, kind=kind(0.0_dp))
        work29 = jets(0, 1, 0, 1)*jets(1, 1, 0, 1)
        work30 = jets(0, 1, 0, 2)*jets(1, 1, 0, 2)
        work31 = work29 + work29 + work30 + work30
        work32 = jets(1, 0, 0, 1)*jets(1, 1, 0, 1)
        work33 = jets(1, 0, 0, 2)*jets(1, 1, 0, 2)
        work34 = jets(0, 1, 0, 1)*jets(2, 0, 0, 1) + jets(0, 1, 0, 2) &
                 *jets(2, 0, 0, 2) + work32 + work33
        work35 = -work34*work12*work14/real(2, kind=kind(0.0_dp))
        work36 = jets(1, 0, 0, 1)*jets(2, 0, 0, 1)
        work37 = jets(1, 0, 0, 2)*jets(2, 0, 0, 2)
        work38 = work36 + work36 + work37 + work37
        work39 = jets(0, 0, 1, 1)*jets(0, 1, 1, 1)
        work40 = jets(0, 0, 1, 2)*jets(0, 1, 1, 2)
        work41 = jets(0, 0, 0, 1)*jets(0, 1, 0, 1)*2 + work39 + &
                 work39 + work40 + work40
        work42 = jets(0, 1, 0, 1)*jets(0, 1, 1, 1)
        work43 = jets(0, 1, 0, 2)*jets(0, 1, 1, 2)
        work44 = jets(0, 0, 1, 1)*jets(0, 2, 0, 1) + jets(0, 0, 1, 2) &
                 *jets(0, 2, 0, 2) + work42 + work43
        work45 = -work44*work14*work15/real(2, kind=kind(0.0_dp))
        work46 = jets(0, 1, 1, 1)*jets(1, 0, 0, 1)
        work47 = jets(0, 1, 1, 2)*jets(1, 0, 0, 2)
        work48 = work19 + work20 + work46 + work47
        work49 = -work48*work12*work15/real(2, kind=kind(0.0_dp))
        work50 = jets(0, 1, 0, 1)*jets(0, 2, 0, 1)
        work51 = jets(0, 1, 0, 2)*jets(0, 2, 0, 2)
        work52 = work50 + work50 + work51 + work51
        work53 = work29 + work30 + jets(0, 2, 0, 1)*jets(1, 0, 0, 1) &
                 + jets(0, 2, 0, 2)*jets(1, 0, 0, 2)
        work54 = -work53*work12*work14/real(2, kind=kind(0.0_dp))
        work55 = work32 + work32 + work33 + work33
        work56 = jets(0, 0, 1, 1)*jets(0, 0, 2, 1)
        work57 = jets(0, 0, 1, 2)*jets(0, 0, 2, 2)
        work58 = jets(0, 0, 0, 1)*jets(0, 0, 1, 1)*2 + work56 + &
                 work56 + work57 + work57
        work59 = work39 + work40 + jets(0, 0, 2, 1)*jets(0, 1, 0, 1) &
                 + jets(0, 0, 2, 2)*jets(0, 1, 0, 2)
        work60 = -work59*work14*work15/real(2, kind=kind(0.0_dp))
        work61 = work16 + work17 + jets(0, 0, 2, 1)*jets(1, 0, 0, 1) &
                 + jets(0, 0, 2, 2)*jets(1, 0, 0, 2)
        work62 = -work61*work12*work15/real(2, kind=kind(0.0_dp))
        work63 = work42 + work42 + work43 + work43
        work64 = work21 + work22 + work46 + work47
        work65 = -work64*work12*work14/real(2, kind=kind(0.0_dp))
        work66 = work25 + work25 + work26 + work26
        energy = work1*work12/real(2, kind=kind(0.0_dp)) + work8 &
                 *work14/real(2, kind=kind(0.0_dp)) + work10*work15/ &
                 real(2, kind=kind(0.0_dp)) + 0
        gh(1) = -aj(1, 2)*work12 - aj(2, 2)*work14 - aj(3, 2)*work15 &
                - work18*work15*work15/real(2, kind=kind(0.0_dp)) - work23 &
                *work14*work15/real(2, kind=kind(0.0_dp)) - work23*work14 &
                *work15/real(2, kind=kind(0.0_dp)) - work27*work12*work15 &
                /real(2, kind=kind(0.0_dp)) - work27*work12*work15/real(2, &
                               kind=kind(0.0_dp)) - work31*work14*work14/real(2, kind= &
                                    kind(0.0_dp)) - work34*work12*work14/real(2, kind= &
                                    kind(0.0_dp)) - work34*work12*work14/real(2, kind= &
                                    kind(0.0_dp)) - work38*work12*work12/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        gh(2) = -aj(1, 3)*work12 - aj(2, 3)*work14 - aj(3, 3)*work15 &
                - work41*work15*work15/real(2, kind=kind(0.0_dp)) - work44 &
                *work14*work15/real(2, kind=kind(0.0_dp)) - work44*work14 &
                *work15/real(2, kind=kind(0.0_dp)) - work48*work12*work15 &
                /real(2, kind=kind(0.0_dp)) - work48*work12*work15/real(2, &
                               kind=kind(0.0_dp)) - work52*work14*work14/real(2, kind= &
                                    kind(0.0_dp)) - work53*work12*work14/real(2, kind= &
                                    kind(0.0_dp)) - work53*work12*work14/real(2, kind= &
                                    kind(0.0_dp)) - work55*work12*work12/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        gh(3) = -aj(1, 4)*work12 - aj(2, 4)*work14 - aj(3, 4)*work15 &
                - work58*work15*work15/real(2, kind=kind(0.0_dp)) - work59 &
                *work14*work15/real(2, kind=kind(0.0_dp)) - work59*work14 &
                *work15/real(2, kind=kind(0.0_dp)) - work61*work12*work15 &
                /real(2, kind=kind(0.0_dp)) - work61*work12*work15/real(2, &
                               kind=kind(0.0_dp)) - work63*work14*work14/real(2, kind= &
                                    kind(0.0_dp)) - work64*work12*work14/real(2, kind= &
                                    kind(0.0_dp)) - work64*work12*work14/real(2, kind= &
                                    kind(0.0_dp)) - work66*work12*work12/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        gh(4) = work12
        gh(5) = work14
        gh(6) = work15
        gj(1, 1) = work3
        gj(2, 1) = work6
        gj(3, 1) = work4
        gj(4, 1) = work5
        gj(5, 1) = work2
        gj(6, 1) = work7
        gj(1, 2) = work38
        gj(2, 2) = work34
        gj(3, 2) = work27
        gj(4, 2) = work31
        gj(5, 2) = work23
        gj(6, 2) = work18
        gj(1, 3) = work55
        gj(2, 3) = work53
        gj(3, 3) = work48
        gj(4, 3) = work52
        gj(5, 3) = work44
        gj(6, 3) = work41
        gj(1, 4) = work66
        gj(2, 4) = work64
        gj(3, 4) = work61
        gj(4, 4) = work63
        gj(5, 4) = work59
        gj(6, 4) = work58

    end subroutine fo_joint_first
    subroutine fo_joint_full(p, aj, jets, gj, energy, gh, hh)
        real(dp), intent(in)::p(3), aj(3, 10), jets(0:3, 0:3, 0:3, 2)
        real(dp), intent(out)::gj(6, 10), energy, gh(6), hh(6, 6)
        real(dp) :: work1, work2, work3, work4, work5, work6, work7, &
                    work8, work9, work10, work11, work12, work13, work14, work15, &
                    work16, work17, work18, work19, work20, work21, work22, &
                    work23, work24, work25, work26, work27, work28, work29, &
                    work30, work31, work32, work33, work34, work35, work36, &
                    work37, work38, work39, work40, work41, work42, work43, &
                    work44, work45, work46, work47, work48, work49, work50, &
                    work51, work52, work53, work54, work55, work56, work57, &
                    work58, work59, work60, work61, work62, work63, work64, &
                    work65, work66, work67, work68, work69, work70, work71, &
                    work72, work73, work74, work75, work76, work77, work78, &
                    work79, work80, work81, work82, work83, work84, work85, &
                    work86, work87, work88, work89, work90, work91, work92, &
                    work93, work94, work95, work96, work97, work98, work99, &
                    work100, work101, work102, work103, work104, work105, work106, &
                    work107, work108, work109, work110, work111, work112, work113, &
                    work114, work115, work116, work117, work118, work119, work120, &
                    work121, work122, work123, work124, work125, work126, work127, &
                    work128, work129, work130, work131, work132, work133, work134, &
                    work135, work136, work137, work138, work139, work140, work141, &
                    work142, work143, work144, work145, work146, work147, work148, &
                    work149, work150, work151, work152, work153, work154, work155, &
                    work156, work157, work158, work159, work160, work161, work162, &
                    work163, work164, work165, work166, work167, work168, work169, &
                    work170, work171, work172, work173, work174, work175, work176, &
                    work177, work178, work179, work180, work181, work182, work183, &
                    work184, work185, work186, work187, work188, work189, work190, &
                    work191, work192, work193, work194, work195, work196, work197, &
                    work198, work199, work200, work201, work202, work203, work204, &
                    work205, work206, work207, work208, work209, work210, work211, &
                    work212, work213, work214, work215, work216, work217, work218, &
                    work219, work220, work221, work222

        work1 = p(1) - aj(1, 1)
        work2 = jets(0, 0, 1, 1)*jets(0, 1, 0, 1) + jets(0, 0, 1, 2)* &
                jets(0, 1, 0, 2)
        work3 = jets(1, 0, 0, 1)*jets(1, 0, 0, 1) + jets(1, 0, 0, 2)* &
                jets(1, 0, 0, 2)
        work4 = jets(0, 0, 1, 1)*jets(1, 0, 0, 1) + jets(0, 0, 1, 2)* &
                jets(1, 0, 0, 2)
        work5 = jets(0, 1, 0, 1)*jets(0, 1, 0, 1) + jets(0, 1, 0, 2)* &
                jets(0, 1, 0, 2)
        work6 = jets(0, 1, 0, 1)*jets(1, 0, 0, 1) + jets(0, 1, 0, 2)* &
                jets(1, 0, 0, 2)
        work7 = jets(0, 0, 0, 1)**2 + jets(0, 0, 1, 1)*jets(0, 0, 1, &
                                                  1) + jets(0, 0, 1, 2)*jets(0, 0, 1, 2)
        work8 = work7*work5 - work2*work2
        work9 = p(2) - aj(2, 1)
        work10 = -work7*work6 + work2*work4
        work11 = p(3) - aj(3, 1)
        work12 = work2*work6 - work4*work5
        work13 = work1*work8/(-work2**2*work3 - work4**2*work5 &
                              - work6**2*work7 + work7*work5*work3 + work2*work4* &
                              work6*2) + work9*work10/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                                     work4*work6*2) + work11*work12/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work14 = work7*work3 - work4*work4
        work15 = -work2*work3 + work4*work6
        work16 = work1*work10/(-work2**2*work3 - work4**2* &
                               work5 - work6**2*work7 + work7*work5*work3 + work2* &
                               work4*work6*2) + work9*work14/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                      work2*work4*work6*2) + work11*work15/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                                        work3 + work2*work4*work6*2) + 0
        work17 = work5*work3 - work6*work6
        work18 = work1*work12/(-work2**2*work3 - work4**2* &
                               work5 - work6**2*work7 + work7*work5*work3 + work2* &
                               work4*work6*2) + work9*work15/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                      work2*work4*work6*2) + work11*work17/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                                        work3 + work2*work4*work6*2) + 0
        work19 = jets(0, 0, 1, 1)*jets(1, 0, 1, 1)
        work20 = jets(0, 0, 1, 2)*jets(1, 0, 1, 2)
        work21 = jets(0, 0, 0, 1)*jets(1, 0, 0, 1)*2 + work19 + &
                 work19 + work20 + work20
        work22 = jets(0, 0, 1, 1)*jets(1, 1, 0, 1)
        work23 = jets(0, 0, 1, 2)*jets(1, 1, 0, 2)
        work24 = jets(0, 1, 0, 1)*jets(1, 0, 1, 1)
        work25 = jets(0, 1, 0, 2)*jets(1, 0, 1, 2)
        work26 = work22 + work23 + work24 + work25
        work27 = -work26*work16*work18/real(2, kind=kind(0.0_dp))
        work28 = jets(1, 0, 0, 1)*jets(1, 0, 1, 1)
        work29 = jets(1, 0, 0, 2)*jets(1, 0, 1, 2)
        work30 = jets(0, 0, 1, 1)*jets(2, 0, 0, 1) + jets(0, 0, 1, 2) &
                 *jets(2, 0, 0, 2) + work28 + work29
        work31 = -work30*work13*work18/real(2, kind=kind(0.0_dp))
        work32 = jets(0, 1, 0, 1)*jets(1, 1, 0, 1)
        work33 = jets(0, 1, 0, 2)*jets(1, 1, 0, 2)
        work34 = work32 + work32 + work33 + work33
        work35 = jets(1, 0, 0, 1)*jets(1, 1, 0, 1)
        work36 = jets(1, 0, 0, 2)*jets(1, 1, 0, 2)
        work37 = jets(0, 1, 0, 1)*jets(2, 0, 0, 1) + jets(0, 1, 0, 2) &
                 *jets(2, 0, 0, 2) + work35 + work36
        work38 = -work37*work13*work16/real(2, kind=kind(0.0_dp))
        work39 = jets(1, 0, 0, 1)*jets(2, 0, 0, 1)
        work40 = jets(1, 0, 0, 2)*jets(2, 0, 0, 2)
        work41 = work39 + work39 + work40 + work40
        work42 = jets(0, 0, 1, 1)*jets(0, 1, 1, 1)
        work43 = jets(0, 0, 1, 2)*jets(0, 1, 1, 2)
        work44 = jets(0, 0, 0, 1)*jets(0, 1, 0, 1)*2 + work42 + &
                 work42 + work43 + work43
        work45 = jets(0, 1, 0, 1)*jets(0, 1, 1, 1)
        work46 = jets(0, 1, 0, 2)*jets(0, 1, 1, 2)
        work47 = jets(0, 0, 1, 1)*jets(0, 2, 0, 1) + jets(0, 0, 1, 2) &
                 *jets(0, 2, 0, 2) + work45 + work46
        work48 = -work47*work16*work18/real(2, kind=kind(0.0_dp))
        work49 = jets(0, 1, 1, 1)*jets(1, 0, 0, 1)
        work50 = jets(0, 1, 1, 2)*jets(1, 0, 0, 2)
        work51 = work22 + work23 + work49 + work50
        work52 = -work51*work13*work18/real(2, kind=kind(0.0_dp))
        work53 = jets(0, 1, 0, 1)*jets(0, 2, 0, 1)
        work54 = jets(0, 1, 0, 2)*jets(0, 2, 0, 2)
        work55 = work53 + work53 + work54 + work54
        work56 = work32 + work33 + jets(0, 2, 0, 1)*jets(1, 0, 0, 1) &
                 + jets(0, 2, 0, 2)*jets(1, 0, 0, 2)
        work57 = -work56*work13*work16/real(2, kind=kind(0.0_dp))
        work58 = work35 + work35 + work36 + work36
        work59 = jets(0, 0, 1, 1)*jets(0, 0, 2, 1)
        work60 = jets(0, 0, 1, 2)*jets(0, 0, 2, 2)
        work61 = jets(0, 0, 0, 1)*jets(0, 0, 1, 1)*2 + work59 + &
                 work59 + work60 + work60
        work62 = work42 + work43 + jets(0, 0, 2, 1)*jets(0, 1, 0, 1) &
                 + jets(0, 0, 2, 2)*jets(0, 1, 0, 2)
        work63 = -work62*work16*work18/real(2, kind=kind(0.0_dp))
        work64 = work19 + work20 + jets(0, 0, 2, 1)*jets(1, 0, 0, 1) &
                 + jets(0, 0, 2, 2)*jets(1, 0, 0, 2)
        work65 = -work64*work13*work18/real(2, kind=kind(0.0_dp))
        work66 = work45 + work45 + work46 + work46
        work67 = work24 + work25 + work49 + work50
        work68 = -work67*work13*work16/real(2, kind=kind(0.0_dp))
        work69 = work28 + work28 + work29 + work29
        work70 = jets(1, 0, 0, 1)*jets(3, 0, 0, 1)
        work71 = jets(1, 0, 0, 2)*jets(3, 0, 0, 2)
        work72 = jets(2, 0, 0, 1)*jets(2, 0, 0, 1)
        work73 = jets(2, 0, 0, 2)*jets(2, 0, 0, 2)
        work74 = work70 + work70 + work71 + work71 + work72 + &
                 work72 + work73 + work73
        work75 = jets(1, 0, 0, 1)*jets(2, 1, 0, 1)
        work76 = jets(1, 0, 0, 2)*jets(2, 1, 0, 2)
        work77 = jets(1, 1, 0, 1)*jets(2, 0, 0, 1)
        work78 = jets(1, 1, 0, 2)*jets(2, 0, 0, 2)
        work79 = jets(0, 1, 0, 1)*jets(3, 0, 0, 1) + jets(0, 1, 0, 2) &
                 *jets(3, 0, 0, 2) + work75 + work76 + work77 + work77 + &
                 work78 + work78
        work80 = jets(1, 0, 0, 1)*jets(2, 0, 1, 1)
        work81 = jets(1, 0, 0, 2)*jets(2, 0, 1, 2)
        work82 = jets(1, 0, 1, 1)*jets(2, 0, 0, 1)
        work83 = jets(1, 0, 1, 2)*jets(2, 0, 0, 2)
        work84 = jets(0, 0, 1, 1)*jets(3, 0, 0, 1) + jets(0, 0, 1, 2) &
                 *jets(3, 0, 0, 2) + work80 + work81 + work82 + work82 + &
                 work83 + work83
        work85 = jets(0, 1, 0, 1)*jets(2, 1, 0, 1)
        work86 = jets(0, 1, 0, 2)*jets(2, 1, 0, 2)
        work87 = jets(1, 1, 0, 1)*jets(1, 1, 0, 1)
        work88 = jets(1, 1, 0, 2)*jets(1, 1, 0, 2)
        work89 = work85 + work85 + work86 + work86 + work87 + &
                 work87 + work88 + work88
        work90 = jets(0, 0, 1, 1)*jets(2, 1, 0, 1)
        work91 = jets(0, 0, 1, 2)*jets(2, 1, 0, 2)
        work92 = jets(0, 1, 0, 1)*jets(2, 0, 1, 1)
        work93 = jets(0, 1, 0, 2)*jets(2, 0, 1, 2)
        work94 = jets(1, 0, 1, 1)*jets(1, 1, 0, 1)
        work95 = jets(1, 0, 1, 2)*jets(1, 1, 0, 2)
        work96 = work90 + work91 + work92 + work93 + work94 + &
                 work94 + work95 + work95
        work97 = jets(0, 0, 1, 1)*jets(2, 0, 1, 1)
        work98 = jets(0, 0, 1, 2)*jets(2, 0, 1, 2)
        work99 = jets(1, 0, 1, 1)*jets(1, 0, 1, 1)
        work100 = jets(1, 0, 1, 2)*jets(1, 0, 1, 2)
        work101 = jets(0, 0, 0, 1)*jets(2, 0, 0, 1)*2 + work97 + &
                  work97 + work98 + work98 + jets(1, 0, 0, 1)*jets(1, 0, 0, 1) &
                  *2 + work99 + work99 + work100 + work100
        work102 = work75 + work75 + work76 + work76 + work77 + &
                  work77 + work78 + work78
        work103 = jets(1, 0, 0, 1)*jets(1, 2, 0, 1)
        work104 = jets(1, 0, 0, 2)*jets(1, 2, 0, 2)
        work105 = work85 + work86 + jets(0, 2, 0, 1)*jets(2, 0, 0, &
                            1) + jets(0, 2, 0, 2)*jets(2, 0, 0, 2) + work103 + work104 &
                  + work87 + work88
        work106 = jets(0, 1, 1, 1)*jets(2, 0, 0, 1)
        work107 = jets(0, 1, 1, 2)*jets(2, 0, 0, 2)
        work108 = jets(1, 0, 0, 1)*jets(1, 1, 1, 1)
        work109 = jets(1, 0, 0, 2)*jets(1, 1, 1, 2)
        work110 = work90 + work91 + work106 + work107 + work108 + &
                  work109 + work94 + work95
        work111 = jets(0, 1, 0, 1)*jets(1, 2, 0, 1)
        work112 = jets(0, 1, 0, 2)*jets(1, 2, 0, 2)
        work113 = jets(0, 2, 0, 1)*jets(1, 1, 0, 1)
        work114 = jets(0, 2, 0, 2)*jets(1, 1, 0, 2)
        work115 = work111 + work111 + work112 + work112 + work113 &
                  + work113 + work114 + work114
        work116 = jets(0, 0, 1, 1)*jets(1, 2, 0, 1)
        work117 = jets(0, 0, 1, 2)*jets(1, 2, 0, 2)
        work118 = jets(0, 1, 0, 1)*jets(1, 1, 1, 1)
        work119 = jets(0, 1, 0, 2)*jets(1, 1, 1, 2)
        work120 = jets(0, 1, 1, 1)*jets(1, 1, 0, 1)
        work121 = jets(0, 1, 1, 2)*jets(1, 1, 0, 2)
        work122 = jets(0, 2, 0, 1)*jets(1, 0, 1, 1)
        work123 = jets(0, 2, 0, 2)*jets(1, 0, 1, 2)
        work124 = work116 + work117 + work118 + work119 + work120 &
                  + work121 + work122 + work123
        work125 = jets(0, 0, 1, 1)*jets(1, 1, 1, 1)
        work126 = jets(0, 0, 1, 2)*jets(1, 1, 1, 2)
        work127 = jets(0, 1, 1, 1)*jets(1, 0, 1, 1)
        work128 = jets(0, 1, 1, 2)*jets(1, 0, 1, 2)
        work129 = jets(0, 0, 0, 1)*jets(1, 1, 0, 1)*2 + work125 + &
                  work125 + work126 + work126 + jets(0, 1, 0, 1)*jets(1, 0, 0, &
                                            1)*2 + work127 + work127 + work128 + work128
        work130 = work103 + work103 + work104 + work104 + work87 &
                  + work87 + work88 + work88
        work131 = work111 + work112 + work113 + work113 + work114 &
                  + work114 + jets(0, 3, 0, 1)*jets(1, 0, 0, 1) + jets(0, 3, 0, &
                                                                     2)*jets(1, 0, 0, 2)
        work132 = jets(0, 2, 1, 1)*jets(1, 0, 0, 1)
        work133 = jets(0, 2, 1, 2)*jets(1, 0, 0, 2)
        work134 = work116 + work117 + work120 + work120 + work121 &
                  + work121 + work132 + work133
        work135 = jets(0, 1, 0, 1)*jets(0, 3, 0, 1)
        work136 = jets(0, 1, 0, 2)*jets(0, 3, 0, 2)
        work137 = jets(0, 2, 0, 1)*jets(0, 2, 0, 1)
        work138 = jets(0, 2, 0, 2)*jets(0, 2, 0, 2)
        work139 = work135 + work135 + work136 + work136 + work137 &
                  + work137 + work138 + work138
        work140 = jets(0, 1, 0, 1)*jets(0, 2, 1, 1)
        work141 = jets(0, 1, 0, 2)*jets(0, 2, 1, 2)
        work142 = jets(0, 1, 1, 1)*jets(0, 2, 0, 1)
        work143 = jets(0, 1, 1, 2)*jets(0, 2, 0, 2)
        work144 = jets(0, 0, 1, 1)*jets(0, 3, 0, 1) + jets(0, 0, 1, 2) &
                  *jets(0, 3, 0, 2) + work140 + work141 + work142 + work142 &
                  + work143 + work143
        work145 = jets(0, 0, 1, 1)*jets(0, 2, 1, 1)
        work146 = jets(0, 0, 1, 2)*jets(0, 2, 1, 2)
        work147 = jets(0, 1, 1, 1)*jets(0, 1, 1, 1)
        work148 = jets(0, 1, 1, 2)*jets(0, 1, 1, 2)
        work149 = jets(0, 0, 0, 1)*jets(0, 2, 0, 1)*2 + work145 + &
                  work145 + work146 + work146 + jets(0, 1, 0, 1)*jets(0, 1, 0, &
                                            1)*2 + work147 + work147 + work148 + work148
        work150 = work80 + work80 + work81 + work81 + work82 + &
                  work82 + work83 + work83
        work151 = work92 + work93 + work106 + work107 + work108 + &
                  work109 + work94 + work95
        work152 = jets(1, 0, 0, 1)*jets(1, 0, 2, 1)
        work153 = jets(1, 0, 0, 2)*jets(1, 0, 2, 2)
        work154 = work97 + work98 + jets(0, 0, 2, 1)*jets(2, 0, 0, &
                            1) + jets(0, 0, 2, 2)*jets(2, 0, 0, 2) + work152 + work153 &
                  + work99 + work100
        work155 = work118 + work118 + work119 + work119 + work120 &
                  + work120 + work121 + work121
        work156 = jets(0, 0, 2, 1)*jets(1, 1, 0, 1)
        work157 = jets(0, 0, 2, 2)*jets(1, 1, 0, 2)
        work158 = jets(0, 1, 0, 1)*jets(1, 0, 2, 1)
        work159 = jets(0, 1, 0, 2)*jets(1, 0, 2, 2)
        work160 = work125 + work126 + work156 + work157 + work158 &
                  + work159 + work127 + work128
        work161 = jets(0, 0, 1, 1)*jets(1, 0, 2, 1)
        work162 = jets(0, 0, 1, 2)*jets(1, 0, 2, 2)
        work163 = jets(0, 0, 2, 1)*jets(1, 0, 1, 1)
        work164 = jets(0, 0, 2, 2)*jets(1, 0, 1, 2)
        work165 = jets(0, 0, 0, 1)*jets(1, 0, 1, 1)*2 + jets(0, 0, &
                              1, 1)*jets(1, 0, 0, 1)*2 + work161 + work161 + work162 + &
                  work162 + work163 + work163 + work164 + work164
        work166 = work108 + work108 + work109 + work109 + work94 &
                  + work94 + work95 + work95
        work167 = work118 + work119 + work120 + work121 + work122 &
                  + work123 + work132 + work133
        work168 = jets(0, 1, 2, 1)*jets(1, 0, 0, 1)
        work169 = jets(0, 1, 2, 2)*jets(1, 0, 0, 2)
        work170 = work125 + work126 + work156 + work157 + work127 &
                  + work128 + work168 + work169
        work171 = work140 + work140 + work141 + work141 + work142 &
                  + work142 + work143 + work143
        work172 = jets(0, 1, 0, 1)*jets(0, 1, 2, 1)
        work173 = jets(0, 1, 0, 2)*jets(0, 1, 2, 2)
        work174 = work145 + work146 + jets(0, 0, 2, 1)*jets(0, 2, 0, &
                            1) + jets(0, 0, 2, 2)*jets(0, 2, 0, 2) + work172 + work173 &
                  + work147 + work148
        work175 = jets(0, 0, 1, 1)*jets(0, 1, 2, 1)
        work176 = jets(0, 0, 1, 2)*jets(0, 1, 2, 2)
        work177 = jets(0, 0, 2, 1)*jets(0, 1, 1, 1)
        work178 = jets(0, 0, 2, 2)*jets(0, 1, 1, 2)
        work179 = jets(0, 0, 0, 1)*jets(0, 1, 1, 1)*2 + jets(0, 0, &
                              1, 1)*jets(0, 1, 0, 1)*2 + work175 + work175 + work176 + &
                  work176 + work177 + work177 + work178 + work178
        work180 = work152 + work152 + work153 + work153 + work99 &
                  + work99 + work100 + work100
        work181 = work158 + work159 + work127 + work127 + work128 &
                  + work128 + work168 + work169
        work182 = work161 + work162 + work163 + work163 + work164 &
                  + work164 + jets(0, 0, 3, 1)*jets(1, 0, 0, 1) + jets(0, 0, 3, &
                                                                     2)*jets(1, 0, 0, 2)
        work183 = work172 + work172 + work173 + work173 + work147 &
                  + work147 + work148 + work148
        work184 = work175 + work176 + work177 + work177 + work178 &
                  + work178 + jets(0, 0, 3, 1)*jets(0, 1, 0, 1) + jets(0, 0, 3, &
                                                                     2)*jets(0, 1, 0, 2)
        work185 = jets(0, 0, 1, 1)*jets(0, 0, 3, 1)
        work186 = jets(0, 0, 1, 2)*jets(0, 0, 3, 2)
        work187 = jets(0, 0, 2, 1)*jets(0, 0, 2, 1)
        work188 = jets(0, 0, 2, 2)*jets(0, 0, 2, 2)
        work189 = jets(0, 0, 0, 1)*jets(0, 0, 2, 1)*2 + jets(0, 0, &
                              1, 1)*jets(0, 0, 1, 1)*2 + work185 + work185 + work186 + &
                  work186 + work187 + work187 + work188 + work188
        work190 = aj(1, 2) + work30*work18 + work37*work16 + &
                  work41*work13
        work191 = aj(2, 2) + work26*work18 + work34*work16 + &
                  work37*work13
        work192 = work190*work191*work10/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work193 = aj(3, 2) + work21*work18 + work26*work16 + &
                  work30*work13
        work194 = work190*work193*work12/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work195 = work191*work193*work15/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work196 = aj(1, 3) + work51*work18 + work56*work16 + &
                  work58*work13
        work197 = aj(2, 3) + work47*work18 + work55*work16 + &
                  work56*work13
        work198 = aj(3, 3) + work44*work18 + work47*work16 + &
                  work51*work13
        work199 = -aj(1, 6)*work13 - aj(2, 6)*work16 - aj(3, 6)* &
                  work18 + work190*work196*work8/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                 work2*work4*work6*2) + work190*work197*work10/(-work2 &
                                  **2*work3 - work4**2*work5 - work6**2*work7 + work7* &
                                 work5*work3 + work2*work4*work6*2) + work190*work198* &
                  work12/(-work2**2*work3 - work4**2*work5 - work6**2* &
                          work7 + work7*work5*work3 + work2*work4*work6*2) + &
                  work196*work191*work10/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                                   work4*work6*2) + work196*work193*work12/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                 work3 + work2*work4*work6*2) + work191*work197*work14 &
                  /(-work2**2*work3 - work4**2*work5 - work6**2*work7 + &
                    work7*work5*work3 + work2*work4*work6*2) + work191* &
                  work198*work15/(-work2**2*work3 - work4**2*work5 - work6 &
                                  **2*work7 + work7*work5*work3 + work2*work4*work6* &
                                2) + work197*work193*work15/(-work2**2*work3 - work4** &
                                 2*work5 - work6**2*work7 + work7*work5*work3 + work2* &
                                   work4*work6*2) + work193*work198*work17/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                work3 + work2*work4*work6*2) - work129*work18*work18*1 &
                  /real(2, kind=kind(0.0_dp)) - work124*work16*work18*2/ &
                  real(2, kind=kind(0.0_dp)) - work110*work13*work18*2/ &
                  real(2, kind=kind(0.0_dp)) - work115*work16*work16*1/ &
                  real(2, kind=kind(0.0_dp)) - work105*work13*work16*2/ &
                  real(2, kind=kind(0.0_dp)) - work102*work13*work13*1/ &
                  real(2, kind=kind(0.0_dp)) + 0
        work200 = aj(1, 4) + work64*work18 + work67*work16 + &
                  work69*work13
        work201 = aj(2, 4) + work62*work18 + work66*work16 + &
                  work67*work13
        work202 = aj(3, 4) + work61*work18 + work62*work16 + &
                  work64*work13
        work203 = -aj(1, 8)*work13 - aj(2, 8)*work16 - aj(3, 8)* &
                  work18 + work190*work200*work8/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                 work2*work4*work6*2) + work190*work201*work10/(-work2 &
                                  **2*work3 - work4**2*work5 - work6**2*work7 + work7* &
                                 work5*work3 + work2*work4*work6*2) + work190*work202* &
                  work12/(-work2**2*work3 - work4**2*work5 - work6**2* &
                          work7 + work7*work5*work3 + work2*work4*work6*2) + &
                  work200*work191*work10/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                                   work4*work6*2) + work200*work193*work12/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                 work3 + work2*work4*work6*2) + work191*work201*work14 &
                  /(-work2**2*work3 - work4**2*work5 - work6**2*work7 + &
                    work7*work5*work3 + work2*work4*work6*2) + work191* &
                  work202*work15/(-work2**2*work3 - work4**2*work5 - work6 &
                                  **2*work7 + work7*work5*work3 + work2*work4*work6* &
                                2) + work201*work193*work15/(-work2**2*work3 - work4** &
                                 2*work5 - work6**2*work7 + work7*work5*work3 + work2* &
                                   work4*work6*2) + work193*work202*work17/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                work3 + work2*work4*work6*2) - work165*work18*work18*1 &
                  /real(2, kind=kind(0.0_dp)) - work160*work16*work18*2/ &
                  real(2, kind=kind(0.0_dp)) - work154*work13*work18*2/ &
                  real(2, kind=kind(0.0_dp)) - work155*work16*work16*1/ &
                  real(2, kind=kind(0.0_dp)) - work151*work13*work16*2/ &
                  real(2, kind=kind(0.0_dp)) - work150*work13*work13*1/ &
                  real(2, kind=kind(0.0_dp)) + 0
        work204 = -work190*work8/(-work2**2*work3 - work4**2* &
                                  work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work191*work10/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work193*work12/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work205 = -work190*work10/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work191*work14/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work193*work15/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work206 = -work190*work12/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work191*work15/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work193*work17/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work207 = work196*work197*work10/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work208 = work196*work198*work12/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work209 = work197*work198*work15/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work210 = -aj(1, 9)*work13 - aj(2, 9)*work16 - aj(3, 9)* &
                  work18 + work196*work200*work8/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                 work2*work4*work6*2) + work196*work201*work10/(-work2 &
                                  **2*work3 - work4**2*work5 - work6**2*work7 + work7* &
                                 work5*work3 + work2*work4*work6*2) + work196*work202* &
                  work12/(-work2**2*work3 - work4**2*work5 - work6**2* &
                          work7 + work7*work5*work3 + work2*work4*work6*2) + &
                  work200*work197*work10/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                                   work4*work6*2) + work200*work198*work12/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                 work3 + work2*work4*work6*2) + work197*work201*work14 &
                  /(-work2**2*work3 - work4**2*work5 - work6**2*work7 + &
                    work7*work5*work3 + work2*work4*work6*2) + work197* &
                  work202*work15/(-work2**2*work3 - work4**2*work5 - work6 &
                                  **2*work7 + work7*work5*work3 + work2*work4*work6* &
                                2) + work201*work198*work15/(-work2**2*work3 - work4** &
                                 2*work5 - work6**2*work7 + work7*work5*work3 + work2* &
                                   work4*work6*2) + work198*work202*work17/(-work2**2* &
                                work3 - work4**2*work5 - work6**2*work7 + work7*work5* &
                                work3 + work2*work4*work6*2) - work179*work18*work18*1 &
                  /real(2, kind=kind(0.0_dp)) - work174*work16*work18*2/ &
                  real(2, kind=kind(0.0_dp)) - work170*work13*work18*2/ &
                  real(2, kind=kind(0.0_dp)) - work171*work16*work16*1/ &
                  real(2, kind=kind(0.0_dp)) - work167*work13*work16*2/ &
                  real(2, kind=kind(0.0_dp)) - work166*work13*work13*1/ &
                  real(2, kind=kind(0.0_dp)) + 0
        work211 = -work196*work8/(-work2**2*work3 - work4**2* &
                                  work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work197*work10/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work198*work12/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work212 = -work196*work10/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work197*work14/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work198*work15/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work213 = -work196*work12/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work197*work15/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work198*work17/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work214 = work200*work201*work10/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work215 = work200*work202*work12/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work216 = work201*work202*work15/(-work2**2*work3 - &
                                   work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                          + work2*work4*work6*2)
        work217 = -work200*work8/(-work2**2*work3 - work4**2* &
                                  work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work201*work10/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work202*work12/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work218 = -work200*work10/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work201*work14/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work202*work15/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work219 = -work200*work12/(-work2**2*work3 - work4**2* &
                                   work5 - work6**2*work7 + work7*work5*work3 + work2* &
                              work4*work6*2) - work201*work15/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                work2*work4*work6*2) - work202*work17/(-work2**2*work3 &
                                 - work4**2*work5 - work6**2*work7 + work7*work5*work3 &
                                                              + work2*work4*work6*2) + 0
        work220 = work10/(-work2**2*work3 - work4**2*work5 - &
                          work6**2*work7 + work7*work5*work3 + work2*work4* &
                          work6*2)
        work221 = work12/(-work2**2*work3 - work4**2*work5 - &
                          work6**2*work7 + work7*work5*work3 + work2*work4* &
                          work6*2)
        work222 = work15/(-work2**2*work3 - work4**2*work5 - &
                          work6**2*work7 + work7*work5*work3 + work2*work4* &
                          work6*2)
        energy = work1*work13/real(2, kind=kind(0.0_dp)) + work9 &
                 *work16/real(2, kind=kind(0.0_dp)) + work11*work18/ &
                 real(2, kind=kind(0.0_dp)) + 0
        gh(1) = -aj(1, 2)*work13 - aj(2, 2)*work16 - aj(3, 2)*work18 &
                - work21*work18*work18/real(2, kind=kind(0.0_dp)) - work26 &
                *work16*work18/real(2, kind=kind(0.0_dp)) - work26*work16 &
                *work18/real(2, kind=kind(0.0_dp)) - work30*work13*work18 &
                /real(2, kind=kind(0.0_dp)) - work30*work13*work18/real(2, &
                               kind=kind(0.0_dp)) - work34*work16*work16/real(2, kind= &
                                    kind(0.0_dp)) - work37*work13*work16/real(2, kind= &
                                    kind(0.0_dp)) - work37*work13*work16/real(2, kind= &
                                    kind(0.0_dp)) - work41*work13*work13/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        gh(2) = -aj(1, 3)*work13 - aj(2, 3)*work16 - aj(3, 3)*work18 &
                - work44*work18*work18/real(2, kind=kind(0.0_dp)) - work47 &
                *work16*work18/real(2, kind=kind(0.0_dp)) - work47*work16 &
                *work18/real(2, kind=kind(0.0_dp)) - work51*work13*work18 &
                /real(2, kind=kind(0.0_dp)) - work51*work13*work18/real(2, &
                               kind=kind(0.0_dp)) - work55*work16*work16/real(2, kind= &
                                    kind(0.0_dp)) - work56*work13*work16/real(2, kind= &
                                    kind(0.0_dp)) - work56*work13*work16/real(2, kind= &
                                    kind(0.0_dp)) - work58*work13*work13/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        gh(3) = -aj(1, 4)*work13 - aj(2, 4)*work16 - aj(3, 4)*work18 &
                - work61*work18*work18/real(2, kind=kind(0.0_dp)) - work62 &
                *work16*work18/real(2, kind=kind(0.0_dp)) - work62*work16 &
                *work18/real(2, kind=kind(0.0_dp)) - work64*work13*work18 &
                /real(2, kind=kind(0.0_dp)) - work64*work13*work18/real(2, &
                               kind=kind(0.0_dp)) - work66*work16*work16/real(2, kind= &
                                    kind(0.0_dp)) - work67*work13*work16/real(2, kind= &
                                    kind(0.0_dp)) - work67*work13*work16/real(2, kind= &
                                    kind(0.0_dp)) - work69*work13*work13/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        gh(4) = work13
        gh(5) = work16
        gh(6) = work18
        gj(1, 1) = work3
        gj(2, 1) = work6
        gj(3, 1) = work4
        gj(4, 1) = work5
        gj(5, 1) = work2
        gj(6, 1) = work7
        gj(1, 2) = work41
        gj(2, 2) = work37
        gj(3, 2) = work30
        gj(4, 2) = work34
        gj(5, 2) = work26
        gj(6, 2) = work21
        gj(1, 3) = work58
        gj(2, 3) = work56
        gj(3, 3) = work51
        gj(4, 3) = work55
        gj(5, 3) = work47
        gj(6, 3) = work44
        gj(1, 4) = work69
        gj(2, 4) = work67
        gj(3, 4) = work64
        gj(4, 4) = work66
        gj(5, 4) = work62
        gj(6, 4) = work61
        gj(1, 5) = work74
        gj(2, 5) = work79
        gj(3, 5) = work84
        gj(4, 5) = work89
        gj(5, 5) = work96
        gj(6, 5) = work101
        gj(1, 6) = work102
        gj(2, 6) = work105
        gj(3, 6) = work110
        gj(4, 6) = work115
        gj(5, 6) = work124
        gj(6, 6) = work129
        gj(1, 7) = work130
        gj(2, 7) = work131
        gj(3, 7) = work134
        gj(4, 7) = work139
        gj(5, 7) = work144
        gj(6, 7) = work149
        gj(1, 8) = work150
        gj(2, 8) = work151
        gj(3, 8) = work154
        gj(4, 8) = work155
        gj(5, 8) = work160
        gj(6, 8) = work165
        gj(1, 9) = work166
        gj(2, 9) = work167
        gj(3, 9) = work170
        gj(4, 9) = work171
        gj(5, 9) = work174
        gj(6, 9) = work179
        gj(1, 10) = work180
        gj(2, 10) = work181
        gj(3, 10) = work182
        gj(4, 10) = work183
        gj(5, 10) = work184
        gj(6, 10) = work189
        hh(1, 1) = -aj(1, 5)*work13 - aj(2, 5)*work16 - aj(3, 5)* &
                   work18 + work190*work190*work8/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                  work2*work4*work6*2) + work192 + work192 + work194 + &
                   work194 + work191*work191*work14/(-work2**2*work3 - &
                                 work4**2*work5 - work6**2*work7 + work7*work5*work3 + &
                                   work2*work4*work6*2) + work195 + work195 + work193* &
                   work193*work17/(-work2**2*work3 - work4**2*work5 - &
                                   work6**2*work7 + work7*work5*work3 + work2*work4* &
                                   work6*2) - work101*work18*work18*1/real(2, kind= &
                                  kind(0.0_dp)) - work96*work16*work18*2/real(2, kind= &
                                  kind(0.0_dp)) - work84*work13*work18*2/real(2, kind= &
                                  kind(0.0_dp)) - work89*work16*work16*1/real(2, kind= &
                                  kind(0.0_dp)) - work79*work13*work16*2/real(2, kind= &
                                  kind(0.0_dp)) - work74*work13*work13*1/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        hh(2, 1) = work199
        hh(3, 1) = work203
        hh(4, 1) = work204
        hh(5, 1) = work205
        hh(6, 1) = work206
        hh(1, 2) = work199
        hh(2, 2) = -aj(1, 7)*work13 - aj(2, 7)*work16 - aj(3, 7)* &
                   work18 + work196*work196*work8/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                  work2*work4*work6*2) + work207 + work207 + work208 + &
                   work208 + work197*work197*work14/(-work2**2*work3 - &
                                 work4**2*work5 - work6**2*work7 + work7*work5*work3 + &
                                   work2*work4*work6*2) + work209 + work209 + work198* &
                   work198*work17/(-work2**2*work3 - work4**2*work5 - &
                                   work6**2*work7 + work7*work5*work3 + work2*work4* &
                                   work6*2) - work149*work18*work18*1/real(2, kind= &
                                 kind(0.0_dp)) - work144*work16*work18*2/real(2, kind= &
                                 kind(0.0_dp)) - work134*work13*work18*2/real(2, kind= &
                                 kind(0.0_dp)) - work139*work16*work16*1/real(2, kind= &
                                 kind(0.0_dp)) - work131*work13*work16*2/real(2, kind= &
                                 kind(0.0_dp)) - work130*work13*work13*1/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        hh(3, 2) = work210
        hh(4, 2) = work211
        hh(5, 2) = work212
        hh(6, 2) = work213
        hh(1, 3) = work203
        hh(2, 3) = work210
        hh(3, 3) = -aj(1, 10)*work13 - aj(2, 10)*work16 - aj(3, 10)* &
                   work18 + work200*work200*work8/(-work2**2*work3 - work4 &
                                      **2*work5 - work6**2*work7 + work7*work5*work3 + &
                                  work2*work4*work6*2) + work214 + work214 + work215 + &
                   work215 + work201*work201*work14/(-work2**2*work3 - &
                                 work4**2*work5 - work6**2*work7 + work7*work5*work3 + &
                                   work2*work4*work6*2) + work216 + work216 + work202* &
                   work202*work17/(-work2**2*work3 - work4**2*work5 - &
                                   work6**2*work7 + work7*work5*work3 + work2*work4* &
                                   work6*2) - work189*work18*work18*1/real(2, kind= &
                                 kind(0.0_dp)) - work184*work16*work18*2/real(2, kind= &
                                 kind(0.0_dp)) - work182*work13*work18*2/real(2, kind= &
                                 kind(0.0_dp)) - work183*work16*work16*1/real(2, kind= &
                                 kind(0.0_dp)) - work181*work13*work16*2/real(2, kind= &
                                 kind(0.0_dp)) - work180*work13*work13*1/real(2, kind= &
                                                                       kind(0.0_dp)) + 0
        hh(4, 3) = work217
        hh(5, 3) = work218
        hh(6, 3) = work219
        hh(1, 4) = work204
        hh(2, 4) = work211
        hh(3, 4) = work217
        hh(4, 4) = work8/(-work2**2*work3 - work4**2*work5 - &
                          work6**2*work7 + work7*work5*work3 + work2*work4* &
                          work6*2)
        hh(5, 4) = work220
        hh(6, 4) = work221
        hh(1, 5) = work205
        hh(2, 5) = work212
        hh(3, 5) = work218
        hh(4, 5) = work220
        hh(5, 5) = work14/(-work2**2*work3 - work4**2*work5 - &
                           work6**2*work7 + work7*work5*work3 + work2*work4* &
                           work6*2)
        hh(6, 5) = work222
        hh(1, 6) = work206
        hh(2, 6) = work213
        hh(3, 6) = work219
        hh(4, 6) = work221
        hh(5, 6) = work222
        hh(6, 6) = work17/(-work2**2*work3 - work4**2*work5 - &
                           work6**2*work7 + work7*work5*work3 + work2*work4* &
                           work6*2)

    end subroutine fo_joint_full
    subroutine fo_joint_flow(energy, e0, tj, gh, hh, kenergy, f, df)
        real(dp), intent(in)::energy, e0, tj(10), gh(6), hh(6, 6)
        real(dp), intent(out)::kenergy, f(6), df(6, 6)
        real(dp) :: work1, work2, work3, work4, work5, work6, work7, &
                    work8, work9, work10, work11, work12, work13, work14, work15, &
                    work16, work17, work18, work19, work20, work21, work22

        work1 = energy - e0
        work2 = gh(4)*tj(2)
        work3 = gh(4)*tj(3)
        work4 = gh(4)*tj(4)
        work5 = gh(5)*tj(2)
        work6 = gh(5)*tj(3)
        work7 = gh(5)*tj(4)
        work8 = gh(6)*tj(2)
        work9 = gh(6)*tj(3)
        work10 = gh(6)*tj(4)
        work11 = gh(1)*tj(2)
        work12 = gh(1)*tj(3)
        work13 = gh(2)*tj(2)
        work14 = tj(6)*work1
        work15 = gh(1)*tj(4)
        work16 = gh(3)*tj(2)
        work17 = tj(8)*work1
        work18 = gh(2)*tj(3)
        work19 = gh(2)*tj(4)
        work20 = gh(3)*tj(3)
        work21 = tj(9)*work1
        work22 = gh(3)*tj(4)
        kenergy = tj(1)*work1
        f(1) = gh(4)*tj(1)
        df(1, 1) = work2 + hh(4, 1)*tj(1)
        df(1, 2) = work3 + hh(4, 2)*tj(1)
        df(1, 3) = work4 + hh(4, 3)*tj(1)
        df(1, 4) = hh(4, 4)*tj(1)
        df(1, 5) = hh(4, 5)*tj(1)
        df(1, 6) = hh(4, 6)*tj(1)
        f(2) = gh(5)*tj(1)
        df(2, 1) = work5 + hh(5, 1)*tj(1)
        df(2, 2) = work6 + hh(5, 2)*tj(1)
        df(2, 3) = work7 + hh(5, 3)*tj(1)
        df(2, 4) = hh(5, 4)*tj(1)
        df(2, 5) = hh(5, 5)*tj(1)
        df(2, 6) = hh(5, 6)*tj(1)
        f(3) = gh(6)*tj(1)
        df(3, 1) = work8 + hh(6, 1)*tj(1)
        df(3, 2) = work9 + hh(6, 2)*tj(1)
        df(3, 3) = work10 + hh(6, 3)*tj(1)
        df(3, 4) = hh(6, 4)*tj(1)
        df(3, 5) = hh(6, 5)*tj(1)
        df(3, 6) = hh(6, 6)*tj(1)
        f(4) = -(gh(1)*tj(1) + tj(2)*work1)
        df(4, 1) = -(work11 + work11 + hh(1, 1)*tj(1) + tj(5)* &
                     work1)
        df(4, 2) = -(work12 + work13 + hh(1, 2)*tj(1) + work14)
        df(4, 3) = -(work15 + work16 + hh(1, 3)*tj(1) + work17)
        df(4, 4) = -(work2 + hh(1, 4)*tj(1))
        df(4, 5) = -(work5 + hh(1, 5)*tj(1))
        df(4, 6) = -(work8 + hh(1, 6)*tj(1))
        f(5) = -(gh(2)*tj(1) + tj(3)*work1)
        df(5, 1) = -(work12 + work13 + hh(2, 1)*tj(1) + work14)
        df(5, 2) = -(work18 + work18 + hh(2, 2)*tj(1) + tj(7)* &
                     work1)
        df(5, 3) = -(work19 + work20 + hh(2, 3)*tj(1) + work21)
        df(5, 4) = -(work3 + hh(2, 4)*tj(1))
        df(5, 5) = -(work6 + hh(2, 5)*tj(1))
        df(5, 6) = -(work9 + hh(2, 6)*tj(1))
        f(6) = -(gh(3)*tj(1) + tj(4)*work1)
        df(6, 1) = -(work15 + work16 + hh(3, 1)*tj(1) + work17)
        df(6, 2) = -(work19 + work20 + hh(3, 2)*tj(1) + work21)
        df(6, 3) = -(work22 + work22 + hh(3, 3)*tj(1) + tj(10)* &
                     work1)
        df(6, 4) = -(work4 + hh(3, 4)*tj(1))
        df(6, 5) = -(work7 + hh(3, 5)*tj(1))
        df(6, 6) = -(work10 + hh(3, 6)*tj(1))

    end subroutine fo_joint_flow
end module fo_joint_kernel
