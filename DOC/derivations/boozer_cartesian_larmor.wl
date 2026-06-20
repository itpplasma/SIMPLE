(* Cartesian Larmor displacement for the CP guiding-center <-> particle map.
   Verifies, symbolically, the three identities the Fortran module
   src/field/boozer_cartesian.f90 relies on:

   (1) The induced metric of the cylindrical embedding equals the pullback form
       boozer_field_metric builds: g_ij = dx/du^i . dx/du^j with
       x = (R cos phi, R sin phi, Z) reproduces gV(3,3)=R^2+R_p^2+Z_p^2 etc.,
       hence g = Jc^T Jc with Jc = dx/du.
   (2) A vector normalized in g has the same Euclidean length under Jc:
       |Jc v|^2 = v^T (Jc^T Jc) v = v^T g v = |v|_g^2.
   (3) Larmor vector normalization and sign: rho = (m/(qc|B|)) (bhat x vperp),
       |rho| = m |vperp|/(qc|B|) = ro0_int |vbar_perp|/|B| with qc=1/ro0_int,
       and the guiding center is x_gc = x_particle - rho. *)

(* --- (1) induced metric = Cartesian Jacobian Gram, for the cyl embedding --- *)
u = {s, th, ph};
R = Rf[s, th, ph]; Z = Zf[s, th, ph]; phi = phf[s, th, ph];
x = {R Cos[phi], R Sin[phi], Z};
Jc = Table[D[x[[a]], u[[k]]], {a, 3}, {k, 3}];
gFromJac = Transpose[Jc].Jc // Simplify;

(* boozer_field_metric form in the geometric angle phi (= varphi_V):
   g_ij = R_i R_j + Z_i Z_j + R^2 phi_i phi_j. *)
Rg = Table[D[R, u[[k]]], {k, 3}];
Zg = Table[D[Z, u[[k]]], {k, 3}];
Pg = Table[D[phi, u[[k]]], {k, 3}];
gPullback = Table[Rg[[i]] Rg[[j]] + Zg[[i]] Zg[[j]] + R^2 Pg[[i]] Pg[[j]],
   {i, 3}, {j, 3}] // Simplify;

Print["(1) g(Jc^T Jc) == cylindrical pullback : ",
  Simplify[gFromJac - gPullback] == ConstantArray[0, {3, 3}]];

(* --- (2) metric norm preserved by Jc (corollary of (1)) --- *)
v = {v1, v2, v3};
lhs = (Jc.v).(Jc.v);                 (* |Jc v|^2 Euclidean *)
rhs = v.(gFromJac.v);                (* v^T g v           *)
Print["(2) |Jc v|^2 == v^T g v      : ", Simplify[lhs - rhs] == 0];

(* --- (3) Larmor normalization and sign --- *)
(* Orthonormal frame (e1,e2,bhat); vperp = vp(e1 cos a + e2 sin a). *)
{e1, e2, bh} = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
vperp = vp (e1 Cos[a] + e2 Sin[a]);
rho = (m/(qc Bmag)) Cross[bh, vperp] // Simplify;
pos = {m > 0, vp > 0, qc > 0, Bmag > 0, ro0int > 0};
Print["(3a) |rho| == m vp/(qc B)    : ",
  Simplify[Sqrt[rho.rho] - m vp/(qc Bmag), pos] == 0];

(* qc = 1/ro0int (e=c=1) gives the SIMPLE gyroradius rho = ro0int vp/B. *)
Print["(3b) |rho| with qc=1/ro0int  : ",
  Simplify[(Sqrt[rho.rho] /. qc -> 1/ro0int) - m vp ro0int/Bmag, pos] == 0];

(* x_gc = x_p - rho is the inverse of x_p = x_gc + rho with the same rho. *)
xp = {xg1, xg2, xg3} + rho;
Print["(3c) (xp - rho) == x_gc      : ",
  Simplify[(xp - rho) - {xg1, xg2, xg3}] == {0, 0, 0}];
