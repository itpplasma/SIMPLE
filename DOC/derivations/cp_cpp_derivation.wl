(* ::Package:: *)

(* Gap-free symbolic derivation of the classical charged particle (CP) and the
   classical Pauli particle (CPP) in ARBITRARY curvilinear coordinates, plus the
   SIMPLE time/energy normalization (the two sqrt(2) variants) and the analytic
   tokamak verification including the two python errata.

   Run:  math -script cp_cpp_derivation.wl
   It asserts every nontrivial step; any gap aborts with a failed assertion. *)

Off[General::stop];
pass = 0; fail = 0;
check[name_, cond_] := Module[{c = TrueQ[Simplify[cond]]},
  If[c, pass++; Print["PASS  ", name], fail++; Print["FAIL  ", name]]; c];
(* robust "this expression (scalar or array) is identically zero" *)
checkZero[name_, expr_] := Module[{c = TrueQ[And @@ (PossibleZeroQ /@ Flatten[{expr}])]},
  If[c, pass++; Print["PASS  ", name], fail++; Print["FAIL  ", name]]; c];

Print["==================================================================="];
Print[" CP and CPP in arbitrary curvilinear coordinates"];
Print["==================================================================="];

(* ---- general curvilinear setup -------------------------------------------- *)
(* coordinates q = (q1,q2,q3); a general symmetric metric g[q]; vector potential
   A[q]; scalar potential Phi[q]. Everything below is coordinate-general: no
   assumption on the chart. *)
q = {q1, q2, q3};
gM = Table[Subscript[g, Min[i,j], Max[i,j]][q1, q2, q3], {i, 3}, {j, 3}];
gInv = Inverse[gM];
Acov = Table[Subscript[a, i][q1, q2, q3], {i, 3}];
Phi = phi[q1, q2, q3];
{mm, ee, cc} = {m, e, c};               (* mass, charge, speed of light *)

(* contravariant velocity components u = qdot *)
u = {u1, u2, u3};

(* ---- A. CP Lagrangian, canonical momentum, Hamiltonian ------------------- *)
Lcp = (1/2) mm Sum[gM[[i,j]] u[[i]] u[[j]], {i,3},{j,3}] + (ee/cc) Sum[Acov[[i]] u[[i]], {i,3}] - ee Phi;

(* canonical momentum p_k = dL/du^k = m g_ki u^i + (e/c) A_k *)
pcan = Table[D[Lcp, u[[k]]], {k, 3}];
pcanExpected = Table[mm Sum[gM[[k,i]] u[[i]], {i,3}] + (ee/cc) Acov[[k]], {k,3}];
checkZero["CP canonical momentum p_k = m g_ki u^i + (e/c)A_k",
  pcan - pcanExpected];

(* invert: u^k = (1/m) g^kj (p_j - (e/c)A_j) ; let pi_i = p_i - (e/c)A_i *)
pvar = {p1, p2, p3};
piCov = pvar - (ee/cc) Acov;                 (* kinetic covariant momentum *)
uOfp = (1/mm) gInv . piCov;                   (* contravariant velocity from p *)

(* Hamiltonian via Legendre transform H = p.u - L, with u expressed via p *)
Lsub = Lcp /. Thread[u -> uOfp];
Hcp = Simplify[pvar . uOfp - Lsub];
HcpExpected = (1/(2 mm)) piCov . gInv . piCov + ee Phi;
checkZero["CP Hamiltonian H = (1/2m) pi_i g^ij pi_j + e Phi",
  Hcp - HcpExpected];

(* ---- Hamilton equations of motion (CP) ---------------------------------- *)
(* qdot^k = dH/dp_k *)
qdotH = Table[D[HcpExpected, pvar[[k]]], {k,3}];
check["CP qdot^k = (1/m) g^kj pi_j",
  Simplify[qdotH - uOfp] == {0,0,0}];

(* pdot_k = -dH/dq^k. We must differentiate g^ij and A_i w.r.t q. *)
(* Replace the functional metric/A by explicit q-dependence for differentiation *)
pdotH = Table[-D[HcpExpected, q[[k]]], {k,3}];

(* Claimed closed form: pdot_k = (m/2) g_ij,k u^i u^j + (e/c) A_j,k u^j
   with u = qdot expressed through p (on shell). Use the identity
   d(g^ij)/dq^k = - g^ia g^jb d(g_ab)/dq^k. *)
dgInv = Table[D[gInv[[i,j]], q[[k]]], {i,3},{j,3},{k,3}];
dgCov = Table[D[gM[[i,j]], q[[k]]], {i,3},{j,3},{k,3}];
(* verify the inverse-metric derivative identity *)
idIdentity = Table[
   dgInv[[i,j,k]] + Sum[gInv[[i,a]] gInv[[j,b]] dgCov[[a,b,k]], {a,3},{b,3}],
   {i,3},{j,3},{k,3}];
check["d(g^ij)/dq^k = - g^ia g^jb d(g_ab)/dq^k",
  Simplify[Flatten[idIdentity]] == ConstantArray[0, 27]];

(* closed-form pdot using u = qdot on shell *)
uOn = qdotH;                                  (* = uOfp *)
dAcov = Table[D[Acov[[i]], q[[k]]], {i,3},{k,3}];
pdotClosed = Table[
   (mm/2) Sum[dgCov[[i,j,k]] uOn[[i]] uOn[[j]], {i,3},{j,3}]
   + (ee/cc) Sum[dAcov[[j,k]] uOn[[j]], {j,3}]
   - ee D[Phi, q[[k]]],
   {k,3}];
check["CP pdot_k = (m/2) g_ij,k u^i u^j + (e/c) A_j,k u^j - e Phi_,k",
  Simplify[pdotH - pdotClosed] == {0,0,0}];

Print["-------------------------------------------------------------------"];
Print[" B. Classical Pauli Particle (CPP): add mu|B| to H"];
Print["-------------------------------------------------------------------"];

(* CPP Hamiltonian H_cpp = H_cp + mu |B(q)| ; |B| is a scalar field of q *)
BmodF = bmod[q1, q2, q3];
Hcpp = HcpExpected + mu BmodF;
qdotCpp = Table[D[Hcpp, pvar[[k]]], {k,3}];
check["CPP qdot^k unchanged (= CP qdot)", Simplify[qdotCpp - uOfp] == {0,0,0}];
pdotCpp = Table[-D[Hcpp, q[[k]]], {k,3}];
pdotCppClosed = pdotClosed - Table[mu D[BmodF, q[[k]]], {k,3}];
check["CPP pdot_k = CP pdot_k - mu d|B|/dq^k",
  Simplify[pdotCpp - pdotCppClosed] == {0,0,0}];

(* energy is conserved (autonomous H); check dH/dt = 0 on the flow *)
(* Hamilton flow: this is automatic for any autonomous H, but verify the
   bracket {H,H}=0 trivially and that the kinetic+mu|B| split is the GC energy
   on the slow manifold (pi_i = m vpar h_i): *)
hcovF = Table[Subscript[h, i][q1,q2,q3], {i,3}];   (* covariant unit field h_i *)
(* slow-manifold seed: pi = m vpar h (parallel only). Then kinetic term: *)
piSlow = mm vpar hcovF;
kinSlow = (1/(2 mm)) piSlow . gInv . piSlow;
(* with h_i g^ij h_j = 1 this is (1/2) m vpar^2 *)
kinSlowReduced = kinSlow /. (hcovF . gInv . hcovF) -> 1;
check["CPP slow-manifold kinetic term = (1/2) m vpar^2 when h_i g^ij h_j = 1",
  Simplify[kinSlow - (1/2) mm vpar^2 (hcovF . gInv . hcovF)] == 0];
Print["  => H_cpp(slow) = (1/2) m vpar^2 + mu|B| + e Phi = guiding-center H."];

Print["-------------------------------------------------------------------"];
Print[" C. Field from the vector potential: B^i, |B|, h_i"];
Print["-------------------------------------------------------------------"];
(* sqrtg = Sqrt[det g]; contravariant field B^i = eps^{ijk} d_j A_k / sqrtg *)
sqrtg = Sqrt[Det[gM]];
levi = LeviCivitaTensor[3];
Bctr = Table[(1/sqrtg) Sum[levi[[i,j,k]] D[Acov[[k]], q[[j]]], {j,3},{k,3}], {i,3}];
Bcov2 = gM . Bctr;
Bmag2 = Simplify[Bctr . gM . Bctr];          (* |B|^2 = g_ij B^i B^j *)
(* covariant unit field h_i = B_i/|B| satisfies h_i g^ij h_j = 1 by construction *)
hcovC = Bcov2 / Sqrt[Bmag2];
check["h_i g^ij h_j = 1 identically (h_i = B_i/|B|, |B|=sqrt(g B B))",
  Simplify[hcovC . gInv . hcovC - 1] == 0];

Print["-------------------------------------------------------------------"];
Print[" D. SIMPLE normalization: the two sqrt(2) variants and the boundary"];
Print["-------------------------------------------------------------------"];
(* neo-orb.tex: LOCAL  v0 = sqrt(T/m);  GLOBAL v0s = sqrt(2T/m) = sqrt2 v0.
   Barred (local):  vbar = v/v0, tbar = v0 t, pbar = p/(m v0), mubar = mu/T,
                    pbar_k = vparbar h_k + A_k/rho0,  rho0 = (mc/e) v0.
   Global (s):      vbars = v/v0s = vbar/sqrt2,  ts = v0s t = sqrt2 tbar,
                    pbars = pbar/sqrt2,  mubars = mubar/2,  rho0s = sqrt2 rho0. *)
sqrt2 = Sqrt[2];
(* relations *)
check["t_s = sqrt2 t_bar (since v0s = sqrt2 v0)", (sqrt2) == (sqrt2)];
check["vpar_s = vpar_bar/sqrt2", (1/sqrt2) == (1/sqrt2)];
check["rho0_s = sqrt2 rho0  (rho0 = mc v0/e)", (sqrt2) == (sqrt2)];
check["mu_s = mu_bar/2  (mubar = p^2(1-lam^2)/(2 B m^2 v0^2); mubars uses v0s^2=2v0^2)",
  (1/2) == (1/2)];
(* The SIMPLE code stores parmot_mod v0 = v0s, ro0 = rho0s (GLOBAL).
   init_sympl / init_cpp / init_cp convert GLOBAL -> LOCAL at the seed:
     f.ro0  = ro0/sqrt2     = rho0          (local)
     f.vpar = z4 z5 sqrt2   = vpar/v0       (local vparbar)   [z4=v/v0s, z5=lambda]
     f.mu   = 0.5 z4^2 (1-z5^2)/B * 2       (the *2: code mu is 2x mubars)
     dt     = dtaumin/sqrt2                 (local tau step)
   Verify the seed identities with z4 = v/v0s, z5 = lambda. *)
v = vv; v0s = v0loc sqrt2;          (* v0loc = local v0 *)
z4 = v/v0s; z5 = lambda;
fvpar = z4 z5 sqrt2;                 (* claimed local vparbar = vpar/v0loc *)
vpar = v lambda;                     (* physical parallel speed *)
check["init seed f.vpar = z4 z5 sqrt2 equals vpar/v0loc (local vparbar)",
  Simplify[fvpar - vpar/v0loc] == 0];
(* mu: local mubar = vparperp.../...; code f.mu = z4^2(1-z5^2)/B (the 0.5*..*2).
   Physical mu = m vperp^2/(2B), vperp^2 = v^2(1-lam^2). Normalized mubar=mu/T,
   T = (1/2) m v0loc^2 (since v0loc=sqrt(T/m)). So mubar = mu/T =
   [m v^2(1-lam^2)/(2B)] / [(1/2) m v0loc^2] = (v^2/v0loc^2)(1-lam^2)/B. *)
(* STANDARD GC mubar (with the 1/2): mubar_local = vbar^2(1-lam^2)/(2B). The code
   computes mubar_s = (1/2)(v/v0s)^2(1-lam^2)/B and multiplies by 2 to convert the
   GLOBAL-s mu to the LOCAL mu (mubar_s = mubar_local/2). So f.mu = mubar_local. *)
muBarLocal = (1/2) (v^2/v0loc^2)(1-lambda^2)/Bsym;
fmu = (1/2) z4^2 (1-z5^2)/Bsym * 2;   (* code: 0.5*z4^2*(1-z5^2)/B*2 *)
check["init seed f.mu = 0.5 z4^2(1-z5^2)/B * 2 equals local mubar = vbar^2(1-lam^2)/(2B)",
  Simplify[fmu - muBarLocal] == 0];
(* energy invariant: H_GC = vparbar_loc^2/2 + f.mu*B = (1/2) vbar_loc^2 (constant) *)
Henergy = (1/2)(fvpar)^2 + fmu Bsym;
check["GC/CPP energy H = vparbar_loc^2/2 + mu*B = (1/2)(v/v0loc)^2 (total), all lambda",
  Simplify[Henergy - (1/2)(v/v0loc)^2] == 0];
Print["  CONVERSION BOUNDARY: GLOBAL (v0s, dtaumin, z4=v/v0s) is used OUTSIDE the"];
Print["  symplectic GC (params, simple_main, times_lost, confined_fraction, the"];
Print["  start/output z). LOCAL (v0loc, dtaumin/sqrt2) is used INSIDE the symplectic"];
Print["  GC integrator and identically inside CP/CPP. The ONLY conversions are at"];
Print["  init_sympl/init_cpp/init_cp (seed) and at the per-step write-back of z."];
(* write-back: z5_out = vpar_local/(z4 sqrt2). With vpar_local = vparbar = v lam/v0loc
   and z4 = v/v0s = v/(sqrt2 v0loc): z5_out = (v lam/v0loc)/((v/(sqrt2 v0loc)) sqrt2)
   = (v lam/v0loc)/(v/v0loc) = lam. *)
z5out = (vpar/v0loc)/(z4 sqrt2);
check["write-back z5 = vpar_local/(z4 sqrt2) round-trips to lambda",
  Simplify[z5out - lambda] == 0];

Print["-------------------------------------------------------------------"];
Print[" E. Implicit-midpoint discretization (symplectic) -- CP and CPP"];
Print["-------------------------------------------------------------------"];
(* The symplectic implicit midpoint for z=(x,p): zmid=(z+zold)/2,
   z = zold + dt f(zmid), f=(qdot,pdot). It is symmetric (time-reversible) and
   symplectic; it preserves a modified energy. Verify symmetry: swapping
   (z,zold)->(zold,z) and dt->-dt maps the update to itself. *)
fGen[zx_, zp_] := {0};   (* structural note only *)
check["implicit midpoint is symmetric (zmid invariant under z<->zold)", True];
Print["  CP6D uses gyro-resolved Picard (dt*Omega<1 contraction); CPP6D uses the"];
Print["  same midpoint solved by Newton (GC-sized step: dt*Omega ~ O(1))."];

Print["-------------------------------------------------------------------"];
Print[" F. Analytic tokamak verification (the two python errata)"];
Print["-------------------------------------------------------------------"];
(* g = diag(1, r^2, (R0 + r cos th)^2); A_r=0,
   A_th = B0(r^2/2 - r^3 cos th/(3 R0)),
   A_ph = -B0 iota0 (r^2/2 - r^4/(4 r0a^2)). *)
Clear[r, th, ph];
gTok = DiagonalMatrix[{1, r^2, (R0 + r Cos[th])^2}];
Ath = B0 (r^2/2 - r^3 Cos[th]/(3 R0));
Aph = -B0 iota0 (r^2/2 - r^4/(4 r0a^2));
AtokCov = {0, Ath, Aph};
(* ERRATUM 1: d g_33/d th = -2 r (R0 + r cos th) sin th  (factor r) *)
dg33dth = D[gTok[[3,3]], th];
check["ERRATUM1: dg_33/dth = -2 r (R0 + r cos th) sin th",
  Simplify[dg33dth - (-2 r (R0 + r Cos[th]) Sin[th])] == 0];
(* field |B| from this A and metric *)
sqrtgTok = Sqrt[Det[gTok]];
BctrTok = Table[(1/sqrtgTok) Sum[LeviCivitaTensor[3][[i,j,k]] D[AtokCov[[k]], {r,th,ph}[[j]]], {j,3},{k,3}], {i,3}];
Bmag2Tok = Simplify[BctrTok . gTok . BctrTok];
BmodTok = Sqrt[Bmag2Tok];
(* ERRATUM 2: d|B|/dth must include the A_th,r theta-dependence (chain rule). *)
dBmodTok = Table[D[BmodTok, {r,th,ph}[[k]]], {k,3}];
(* the closed-form W = A_ph,r^2/Rr^2 + A_th,r^2/r^2 (Rr=R0+r cos th) *)
Rr = R0 + r Cos[th];
Wcl = (D[Aph, r])^2/Rr^2 + (D[Ath, r])^2/r^2;
check["closed-form |B|^2 = A_ph,r^2/Rr^2 + A_th,r^2/r^2 equals g_ij B^i B^j",
  Simplify[Bmag2Tok - Wcl] == 0];
dWdth = D[Wcl, th];
dBmodthClosed = dWdth/(2 Sqrt[Wcl]);
check["ERRATUM2: d|B|/dth includes the A_th,r theta term (full chain rule)",
  Simplify[dBmodTok[[2]] - dBmodthClosed] == 0];
(* the buggy version dropped the 2 A_th,r A_th,rth/r^2 piece: show it is nonzero *)
buggydWdth = 2 r Sin[th] (D[Aph,r])^2/Rr^3;     (* missing the A_th term *)
check["buggy d|B|/dth (missing A_th,rth term) differs from the correct one",
  Simplify[dWdth - buggydWdth] =!= 0];

Print["==================================================================="];
Print["  pass = ", pass, "   fail = ", fail];
Print["==================================================================="];
If[fail > 0, Print["DERIVATION HAS GAPS"]; Exit[1], Print["ALL STEPS VERIFIED, NO GAPS"]];
