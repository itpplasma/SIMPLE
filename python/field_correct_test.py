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

        self.ctr_Br  = 1/J *(self.co_dAph[1] - self.co_dAth[2])
        self.ctr_Bth = 1/J *(-self.co_dAph[0]) 
        self.ctr_Bph = 1/J *(self.co_dAth[0])

        self.B = sqrt(self.ctr_Br**2 + r**2 *self.ctr_Bth**2 + (R0+r*cth)**2 *self.ctr_Bph**2)
        self.dB = array(((r*self.ctr_Bth**2 + r**2*self.ctr_Bth*(-1/J*self.co_d2Aph[0] + 1/J**2*(R0+2*r*cth)*self.co_dAph[0]) + (R0+r*cth)*cth*self.ctr_Bph**2 + (R0+r*cth)**2*self.ctr_Bph*(1/J*self.co_d2Ath[0] - 1/J**2*(R0+2*r*cth)*self.co_dAth[0]))/self.B,
                         (-(R0+r*cth)*r*sth*self.ctr_Bph**2 + (R0+r*cth)**2*self.ctr_Bph*(1/J*self.co_d2Ath[3] + 1/J**2*r**2*sth*self.co_dAth[0]))/self.B,
                         zer))

        self.B1 = B0*(iota0*(r-r**3/a**2)/(R0+r*cth) + (1.0 - r/R0*cth))

        self.co_hr  = self.ctr_Br/self.B
        self.co_hth = r**2 *self.ctr_Bth/self.B
        self.co_hph = (R0+r*cth)**2 *self.ctr_Bph/self.B 



        self.Phie = zer
        self.dPhie = array((zer, zer, zer))
        self.d2Phie = array((zer, zer, zer, zer, zer, zer))


        #self.B1    = B0*(iota0*(r-r**3/a**2)/(R0+r*cth) + (1.0 - r/R0*cth))
        #self.dB   = array((B0*(iota0*(1-3*r**2/a**2)/(R0+r*cth) - iota0*(r-r**3/a**2)/(R0+r*cth)**2 *cth - cth/R0), 
        #                   B0*(iota0*(r-r**3/a**2)/(R0+r*cth)**2 *r*sth + r*sth/R0),
        #                   zer))
        #self.d2B  = array((B0*(-iota0*6*r/a**2/(R0+r*cth) - iota0*(1-3*r**2/a**2)/(R0+r*cth)**2*cth ................... )))
        #self.dB   = array((-B0/R0*cth, B0*r/R0*sth, zer))
        #self.d2B  = array((zer, B0*r/R0*cth, zer, B0/R0*sth, zer, zer))
