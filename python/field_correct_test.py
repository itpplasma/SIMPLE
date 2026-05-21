from numpy import sin, cos, zeros_like, array, sqrt

B0 = 1.0    # magnetic field modulus normalization
iota0 = 1.0 # constant part of rotational transform
a = 0.5     # (equivalent) minor radius
R0 = 1.0    # (equivalent) major radius

class field:
    """ Model tokamak field with exact B = B0*(iota0*(r-r^3/a^2)/(R0+r*cos(th)) + (1-r/R0*cos(th)) """
    def evaluate(self, r, th, ph):
        cth = cos(th)
        sth = sin(th)
        zer = zeros_like(r)
        J = r*(R0+r*cth)

        self.co_Ath = B0*(r**2/2.0 - r**3/(3.0*R0)*cth)
        self.co_dAth = array((B0*(r - r**2/R0*cth), B0*r**3*sth/(3.0*R0), zer))
        self.co_d2Ath = array((B0*(1.0 - 2.0*r/R0*cth),
                            B0*r**3*cth/(3.0*R0),
                            zer,
                            B0*r**2/R0*sth,
                            zer,
                            zer))

        self.co_Aph  = -B0*iota0*(r**2/2.0-r**4/(4.0*a**2))
        self.co_dAph = array((-B0*iota0*(r-r**3/a**2), zer, zer))
        self.co_d2Aph = array((-B0*iota0*(1.0-3.0*r**2/a**2), zer, zer, zer, zer, zer))

        # Scalar potential
        self.Phie = zer
        self.dPhie = array((zer, zer, zer))
        self.d2Phie = array((zer, zer, zer, zer, zer, zer))


        # Script by ChatGPT, verified by hand
        # ------------------------------------------------------------------
        # Magnetic-field modulus and unit-field covariant components
        # ------------------------------------------------------------------

        # Short notation
        D = R0 + r*cth          # major-radius factor
        J = r*D

        # Contravariant magnetic field components
        self.ctr_Br  = zer
        self.ctr_Bth = -self.co_dAph[0] / J
        self.ctr_Bph =  self.co_dAth[0] / J

        # It is useful to write the physical B-components as
        #
        #   U = r B^theta
        #   V = D B^phi
        #
        # Then |B|^2 = U^2 + V^2.
        #
        # Since
        #   B^theta = S / J,
        #   B^phi   = Q / J,
        #
        # we have
        #   U = S / D,
        #   V = Q / r.

        S   = B0*iota0*(r - r**3/a**2)
        Sr  = B0*iota0*(1.0 - 3.0*r**2/a**2)
        Srr = B0*iota0*(-6.0*r/a**2)

        Q   = B0*(r - r**2/R0*cth)
        Qr  = B0*(1.0 - 2.0*r/R0*cth)
        Qrr = B0*(-2.0/R0*cth)

        Qt  = B0*(r**2/R0*sth)
        Qtt = B0*(r**2/R0*cth)
        Qrt = B0*(2.0*r/R0*sth)

        Dr  = cth
        Dt  = -r*sth
        Drr = zer
        Dtt = -r*cth
        Drt = -sth

        def div2(A, Ar, At, Arr, Att, Art,
                 C, Cr, Ct, Crr, Ctt, Crt):
            """
            Derivatives up to second order of F = A / C.

            Returns:
                F, Fr, Ft, Frr, Ftt, Frt
            """
            F = A / C

            Fr = Ar/C - A*Cr/C**2
            Ft = At/C - A*Ct/C**2

            Frr = Arr/C - (2.0*Ar*Cr + A*Crr)/C**2 + 2.0*A*Cr**2/C**3
            Ftt = Att/C - (2.0*At*Ct + A*Ctt)/C**2 + 2.0*A*Ct**2/C**3

            Frt = (
                Art/C
                - (Ar*Ct + At*Cr + A*Crt)/C**2
                + 2.0*A*Cr*Ct/C**3
            )

            return F, Fr, Ft, Frr, Ftt, Frt

        # U = S / D
        U, Ur, Ut, Urr, Utt, Urt = div2(
            S, Sr, zer, Srr, zer, zer,
            D, Dr, Dt, Drr, Dtt, Drt
        )

        # V = Q / r
        V, Vr, Vt, Vrr, Vtt, Vrt = div2(
            Q, Qr, Qt, Qrr, Qtt, Qrt,
            r, 1.0, zer, zer, zer, zer
        )

        # Magnetic-field modulus B = sqrt(U^2 + V^2)
        G = U**2 + V**2

        Gr = 2.0*(U*Ur + V*Vr)
        Gt = 2.0*(U*Ut + V*Vt)

        Grr = 2.0*(Ur**2 + U*Urr + Vr**2 + V*Vrr)
        Gtt = 2.0*(Ut**2 + U*Utt + Vt**2 + V*Vtt)
        Grt = 2.0*(Ur*Ut + U*Urt + Vr*Vt + V*Vrt)

        self.B = sqrt(G)

        Br  = Gr/(2.0*self.B)
        Bt  = Gt/(2.0*self.B)

        Brr = Grr/(2.0*self.B) - Gr**2/(4.0*self.B**3)
        Btt = Gtt/(2.0*self.B) - Gt**2/(4.0*self.B**3)
        Brt = Grt/(2.0*self.B) - Gr*Gt/(4.0*self.B**3)

        self.dB = array((Br, Bt, zer))
        self.d2B = array((Brr, Btt, zer, Brt, zer, zer))

        # ------------------------------------------------------------------
        # Covariant components h_i = g_ij B^j / |B|
        # ------------------------------------------------------------------

        self.co_hr = zer
        self.co_dhr = array((zer, zer, zer))
        self.co_d2hr = array((zer, zer, zer, zer, zer, zer))

        # h_theta numerator:
        #
        #   h_theta = r^2 B^theta / B
        #           = r U / B
        #
        Nth   = r*U
        Nthr  = U + r*Ur
        Ntht  = r*Ut
        Nthrr = 2.0*Ur + r*Urr
        Nthtt = r*Utt
        Nthrt = Ut + r*Urt

        hth, hthr, htht, hthrr, hthtt, hthrt = div2(
            Nth, Nthr, Ntht, Nthrr, Nthtt, Nthrt,
            self.B, Br, Bt, Brr, Btt, Brt
        )

        self.co_hth = hth
        self.co_dhth = array((hthr, htht, zer))
        self.co_d2hth = array((hthrr, hthtt, zer, hthrt, zer, zer))

        # h_phi numerator:
        #
        #   h_phi = D^2 B^phi / B
        #         = D V / B
        #
        Nph = D*V

        Nphr = Dr*V + D*Vr
        Npht = Dt*V + D*Vt

        Nphrr = Drr*V + 2.0*Dr*Vr + D*Vrr
        Nphtt = Dtt*V + 2.0*Dt*Vt + D*Vtt
        Nphrt = Drt*V + Dr*Vt + Dt*Vr + D*Vrt

        hph, hphr, hpht, hphrr, hphtt, hphrt = div2(
            Nph, Nphr, Npht, Nphrr, Nphtt, Nphrt,
            self.B, Br, Bt, Brr, Btt, Brt
        )

        self.co_hph = hph
        self.co_dhph = array((hphr, hpht, zer))
        self.co_d2hph = array((hphrr, hphtt, zer, hphrt, zer, zer))





'''

        self.co_Ath = B0*(r**2/2.0 - r**3/(3.0*R0)*cth)
        self.co_dAth = array((B0*(r - r**2/R0*cth), B0*r**3*sth/(3.0*R0), zer))
        self.co_d2Ath = array((B0*(1.0 - 2.0*r/R0*cth),
                            B0*r**3*cth/(3.0*R0),
                            zer,
                            B0*r**2/R0*sth,
                            zer,
                            zer))

        self.co_Aph  = -B0*iota0*(r**2/2.0-r**4/(4.0*a**2))
        self.co_dAph = array((-B0*iota0*(r-r**3/a**2), zer, zer))
        self.co_d2Aph = array((-B0*iota0*(1.0-3.0*r**2/a**2), zer, zer, zer, zer, zer))

        self.ctr_Br  = 1/J *(self.co_dAph[1] - self.co_dAth[2]) # zero
        self.ctr_Bth = 1/J *(-self.co_dAph[0])
        self.ctr_Bph = 1/J *(self.co_dAth[0])

        self.B = sqrt(self.ctr_Br**2 + r**2 *self.ctr_Bth**2 + (R0+r*cth)**2 *self.ctr_Bph**2)
        self.dB = array(((r*self.ctr_Bth**2 + r**2*self.ctr_Bth*(-1/J*self.co_d2Aph[0] + 1/J**2*(R0+2*r*cth)*self.co_dAph[0]) + (R0+r*cth)*cth*self.ctr_Bph**2 + (R0+r*cth)**2*self.ctr_Bph*(1/J*self.co_d2Ath[0] - 1/J**2*(R0+2*r*cth)*self.co_dAth[0]))/self.B,
                         (-r**2*self.ctr_Bth*(1/J**2 *r**2*sth*self.co_dAph[0] + 1/J*self.co_d2Aph[3]) - (R0+r*cth)*r*sth*self.ctr_Bph**2 + (R0+r*cth)**2*self.ctr_Bph*(1/J*self.co_d2Ath[3] + 1/J**2*r**2*sth*self.co_dAth[0]))/self.B,
                         zer))

        self.B1 = B0*(iota0*(r-r**3/a**2)/(R0+r*cth) + (1.0 - r/R0*cth))

        self.co_hr  = self.ctr_Br/self.B # zero
        self.co_dhr = array((zer, zer, zer))
        self.co_d2hr = array((zer, zer, zer, zer, zer, zer))

        self.co_hth = r**2 *self.ctr_Bth/self.B

        self.co_dhth = array((2*r*self.ctr_Bth/self.B + r**2*(1/J**2 *r**2*sth*self.co_dAph[0] + 1/J*self.co_d2Aph[0]) - r**2*self.ctr_Bth*(1/J**2 *r**2*sth*self.co_dAph[0] + 1/J*self.co_d2Aph[0]) - r**2*self.ctr_Bth*(-r*self.ctr_Bth**2 - (R0+r*cth)*cth*self.ctr_Bph**2 - (R0+r*cth)**2*self.ctr_Bph*(1/J*self.co_d2Ath[0] - 1/J**2*(R0+2*r*cth)*self.co_dAth[0]))/self.B),
                            zer,
                            zer))
        self.co_hph = (R0+r*cth)**2 *self.ctr_Bph/self.B



        self.Phie = zer
        self.dPhie = array((zer, zer, zer))
        self.d2Phie = array((zer, zer, zer, zer, zer, zer))
'''
