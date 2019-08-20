"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""
from numpy import sin, cos, zeros_like, array

B0 = 1.0    # magnetic field modulus normalization
iota0 = 1.0 # constant part of rotational transform
a = 0.5     # (equivalent) minor radius
R0 = 1.0    # (equivalent) major radius

class field:
    """ Model tokamak field with B ~ B0*(1-r/R0*cos(th)) """
    def evaluate(self, r, th, ph):
        cth = cos(th)
        sth = sin(th)
        zer = zeros_like(r)

        self.Ath = B0*(r**2/2.0 - r**3/(3.0*R0)*cth)
        self.dAth = array((B0*(r - r**2/R0*cth), B0*r**3*sth/(3.0*R0), zer))
        self.d2Ath = array((B0*(1.0 - 2.0*r/R0*cth), 
                            B0*r**3*cth/(3.0*R0), 
                            zer,
                            B0*r**2/R0*sth, 
                            zer, 
                            zer))

        self.Aph  = -B0*iota0*(r**2/2.0-r**4/(4.0*a**2))
        self.dAph = array((-B0*iota0*(r-r**3/a**2), zer, zer))
        self.d2Aph = array((-B0*iota0*(1.0-3.0*r**2/a**2), zer, zer, zer, zer, zer))

        self.hth  = iota0*(1.0-r**2/a**2)*r**2/R0
        self.dhth = array(((2.0*iota0*r*(a**2-2.0*r**2))/(a**2*R0), zer, zer))
        self.d2hth = array(((2.0*iota0*(a**2-6.0*r**2))/(a**2*R0), zer, zer, zer, zer, zer))

        self.hph  = R0 + r*cth
        self.dhph = array((cth, -r*sth, zer))
        self.d2hph = array((zer, -r*cth, zer, -sth, zer) )

        self.B    = B0*(1.0 - r/R0*cth)
        self.dB   = array((-B0/R0*cth, B0*r/R0*sth, zer))
        self.d2B  = array((zer, B0*r/R0*cth, zer, B0/R0*sth, zer, zer))

        self.Phie = zer
        self.dPhie = array((zer, zer, zer))
        self.d2Phie = array((zer, zer, zer, zer, zer, zer))
