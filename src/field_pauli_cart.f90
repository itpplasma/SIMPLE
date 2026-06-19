module field_pauli_cart
  ! Analytic Cartesian circular-tokamak field for the genuine 6D classical Pauli
  ! particle (orbit_cpp_pauli). The vector potential is exact and the field is
  ! divergence-free by construction (B = curl A), so the canonical Hamiltonian
  !   H = |p - (q/c) A|^2/(2m) + mu |B|
  ! is well posed. This is the SHARED equilibrium of field_can_test and
  ! orbit_full_tokamak (R0, a, B0, iota0); near the axis B_phi = B0 R0/R and the
  ! poloidal field grows like iota0 r/R0, matching field_can_test to leading
  ! order. Unlike orbit_full_tokamak's near-axis ansatz, this B is exactly
  ! solenoidal because it is defined through A.
  !
  !   psi   = B0 iota0/(2 R0) ((R-R0)^2 + Z^2),  A_phi = psi/R
  !   A_z   = -B0 R0/2 ln(R^2),  A_R = 0  (R = sqrt(x^2+y^2))
  !
  ! Everything needed by an implicit-symplectic step with an ANALYTIC Jacobian is
  ! returned in one pass: A, grad A, Hess A, |B|, grad|B|, Hess|B|. The routine
  ! is pure and !$acc routine seq, no procedure pointers or class() dispatch, so
  ! it inlines into the device kernel.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: pauli_field_params_t, eval_pauli_field_cart

  ! Second-derivative packing index m over symmetric pairs (j,k):
  !   1:(x,x) 2:(x,y) 3:(x,z) 4:(y,y) 5:(y,z) 6:(z,z)
  type :: pauli_field_params_t
    real(dp) :: B0    = 1.0_dp
    real(dp) :: R0    = 1.0_dp
    real(dp) :: iota0 = 1.0_dp
    real(dp) :: a     = 0.5_dp   ! minor radius, carried for slow-manifold setup
  end type pauli_field_params_t

contains

  ! Evaluate A, grad A, Hess A, |B|, grad|B|, Hess|B| at Cartesian xv = (x,y,z).
  ! dA(i,j)   = dA_i/dx_j ;  d2A(i,m)  = d2A_i/(dx_j dx_k) for pair m.
  ! dBmod(j)  = d|B|/dx_j  ; d2Bmod(m) = d2|B|/(dx_j dx_k) for pair m.
  pure subroutine eval_pauli_field_cart(p, xv, Avec, dA, d2A, Bvec, Bmod, &
      dBmod, d2Bmod)
    !$acc routine seq
    type(pauli_field_params_t), intent(in) :: p
    real(dp), intent(in)  :: xv(3)
    real(dp), intent(out) :: Avec(3), dA(3,3), d2A(3,6)
    real(dp), intent(out) :: Bvec(3), Bmod, dBmod(3), d2Bmod(6)
    real(dp) :: x, y, z, R0, B0, iota0
    real(dp) :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, &
      t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, &
      t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, &
      t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, &
      t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, &
      t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95, &
      t96, t97, t98, t99, t100, t101, t102, t103, t104, t105, t106, t107, t108, t109, &
      t110, t111, t112, t113, t114, t115, t116, t117, t118, t119, t120, t121, t122, t123, &
      t124, t125, t126, t127, t128

    x = xv(1); y = xv(2); z = xv(3)
    R0 = p%R0; B0 = p%B0; iota0 = p%iota0

    t0 = x**2
    t1 = y**2
    t2 = t0 + t1
    t3 = 1d0/t2
    t4 = z**2
    t5 = sqrt(t2)
    t6 = R0 - t5
    t7 = t4 + t6**2
    t8 = (1.0d0/2.0d0)*t3*t7
    t9 = iota0/R0
    t10 = B0*t9
    t11 = t10*y
    t12 = t10*x
    t13 = t6/t2**(3.0d0/2.0d0)
    t14 = t2**(-2)
    t15 = t14*t7
    t16 = t11*x*(t13 + t15)
    t17 = -t8
    t18 = t1*t13
    t19 = t1*t15 + t18
    t20 = B0*t3
    t21 = t9*z
    t22 = t21*y
    t23 = t0*t13
    t24 = t0*t15 + t23
    t25 = t21*x
    t26 = R0*x
    t27 = R0*y
    t28 = t6/t2**(5.0d0/2.0d0)
    t29 = 4*t28
    t30 = t0*t3
    t31 = 4*t30
    t32 = -t6/t5
    t33 = t0*t29 + t15*(t31 - 1) + t3*(t23 + t30 + t32)
    t34 = -t13 - t15
    t35 = t1*t14
    t36 = t2**(-3)
    t37 = t36*t7
    t38 = 4*t37
    t39 = t1*t38
    t40 = 5*t28
    t41 = t1*t40
    t42 = t35 + t39 + t41
    t43 = 2*t14
    t44 = B0*t43
    t45 = t25*t44*y
    t46 = t1*t3
    t47 = 4*t46
    t48 = t1*t29 + t15*(t47 - 1) + t3*(t18 + t32 + t46)
    t49 = -2*t13 - 2*t15
    t50 = 2*t46
    t51 = t50 - 1
    t52 = t20*t21
    t53 = t20*t9
    t54 = t0*t14
    t55 = t0*t38
    t56 = t0*t40
    t57 = t54 + t55 + t56
    t58 = 2*t30
    t59 = t58 - 1
    t60 = R0*t20
    t61 = t25 + t27
    t62 = -t22
    t63 = t26 + t62
    t64 = t19 + t24 - t3*t7
    t65 = -t64
    t66 = iota0**2/R0**2
    t67 = t14*t61**2 + t14*t63**2 + t64**2*t66
    t68 = sqrt(t67)
    t69 = t3*y
    t70 = 2*t25*t69
    t71 = -R0*t58 + R0 + t70
    t72 = t3*x
    t73 = 2*t27*t72
    t74 = t21*t58 - t21 + t73
    t75 = t14*t61
    t76 = 4*t15
    t77 = 4*t13
    t78 = t42 + t57 - t76 - t77
    t79 = t64*t66
    t80 = t78*t79
    t81 = -t14*t63*t71 + t74*t75 + t80*x
    t82 = 1d0/t68
    t83 = B0*t82
    t84 = -t21*t50 + t21 + t73
    t85 = R0*t50 - R0 + t70
    t86 = t14*t63*t84 + t75*t85
    t87 = t61*t72
    t88 = t63*t69
    t89 = t30 + t46
    t90 = t89 - 1
    t91 = 2*t90
    t92 = t53*t82
    t93 = x**3
    t94 = 4*t3
    t95 = R0*t94
    t96 = t22*t31
    t97 = 2*t36
    t98 = t63*t97
    t99 = t21*t94
    t100 = t27*t31 - t27
    t101 = t61*t97
    t102 = t66*t78**2
    t103 = x**4
    t104 = 9*t36
    t105 = 35*t28
    t106 = 28*t37
    t107 = t7/t2**4
    t108 = 24*t107
    t109 = t6/t2**(7.0d0/2.0d0)
    t110 = 33*t109
    t111 = t0*t1
    t112 = t104*t111 + t108*t111 + t110*t111 + t76 + t77
    t113 = 1d0/t67
    t114 = t25*t47
    t115 = t26*t47 - t26
    t116 = x*y
    t117 = 3*t36
    t118 = 11*t109
    t119 = 8*t107
    t120 = t80*y + t86
    t121 = t113*t81
    t122 = t116*t43
    t123 = 8*t21*t64*t90
    t124 = t78*t91
    t125 = t64*t91
    t126 = t125*t21 + t87 - t88
    t127 = y**3
    t128 = y**4
    Avec(1) = -t11*t8
    Avec(2) = t12*t8
    Avec(3) = -1.0d0/2.0d0*B0*R0*log(t2)
    dA(1,1) = t16
    dA(1,2) = t10*(t17 + t19)
    dA(1,3) = -t20*t22
    dA(2,1) = t10*(-t17 - t24)
    dA(2,2) = -t16
    dA(2,3) = t20*t25
    dA(3,1) = -t20*t26
    dA(3,2) = -t20*t27
    dA(3,3) = 0
    d2A(1,1) = -t11*t33
    d2A(1,2) = -t12*(t34 + t42)
    d2A(1,3) = t45
    d2A(1,4) = -t11*(t48 + t49)
    d2A(1,5) = t51*t52
    d2A(1,6) = -t53*y
    d2A(2,1) = t12*(t33 + t49)
    d2A(2,2) = t11*(t34 + t57)
    d2A(2,3) = -t52*t59
    d2A(2,4) = t12*t48
    d2A(2,5) = -t45
    d2A(2,6) = t53*x
    d2A(3,1) = t59*t60
    d2A(3,2) = t27*t44*x
    d2A(3,3) = 0
    d2A(3,4) = t51*t60
    d2A(3,5) = 0
    d2A(3,6) = 0
    Bvec(1) = -t20*t61
    Bvec(2) = t20*t63
    Bvec(3) = t10*t65
    Bmod = B0*t68
    dBmod(1) = -t81*t83
    dBmod(2) = -t83*(-t65*t66*t78*y + t86)
    dBmod(3) = -t92*(t21*t65*t91 - t87 + t88)
    d2Bmod(1) = t83*(t0*t102 + t101*(t100 - 3*t25 + t93*t99) - t113*t81**2 + t14*t71**2 &
      + t14*t74**2 + t79*(-t0*t105 - t0*t106 + t103*t104 + t103*t108 + &
      t103*t110 + t112 - t35 - t39 - t41 - 7*t54) + t98*(t22 - 3*t26 + &
      t93*t95 - t96))
    d2Bmod(2) = t83*(t101*(t115 + t62 + t96) + t102*t116 + 3*t116*t79*(t0*t117 + t0*t118 &
      + t0*t119 + t1*t117 + t1*t118 + t1*t119 - 10*t28 - 8*t37 - t43) - &
      t120*t121 - t14*t71*t84 + t14*t74*t85 + t98*(t100 - t114 + t25))
    d2Bmod(3) = -t92*(-t121*t126 - t122*t63 + t123*t72 + t124*t25 + t3*t59*t61 + t69*t71 &
      + t72*t74)
    d2Bmod(4) = t83*(t1*t102 + t101*(t114 + t127*t95 - t25 - 3*t27) - t113*t120**2 + t14 &
      *t84**2 + t14*t85**2 + t79*(-t1*t105 - t1*t106 + t104*t128 + t108 &
      *t128 + t110*t128 + t112 - 7*t35 - t54 - t55 - t56) + t98*(t115 - &
      t127*t99 + 3*t22))
    d2Bmod(5) = -t92*(-t113*t120*t126 + t122*t61 + t123*t69 + t124*t22 - t3*t51*t63 - &
      t69*t84 + t72*t85)
    d2Bmod(6) = t20*t66*t82*(-t113*t126**2*t3 + t125 + t4*t90**2*t94 + t89)
  end subroutine eval_pauli_field_cart

end module field_pauli_cart
