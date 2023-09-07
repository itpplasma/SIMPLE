"""
Created: Tue Aug 20 15:07:28 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""

import numpy as np
from fffi import fortran_library, fortran_module

libneo_orb = fortran_library('neo_orb', path='/u/calbert/build/NEO-ORB',
                             compiler={'name': 'ifort', 'version': '18.0.3'})

neoorb = fortran_module(libneo_orb, 'neo_orb_global')
neoorb.fdef("""
  subroutine init_field(ans_s, ans_tp, amultharm, aintegmode)
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
  end
""")

libneo_orb.compile()
neoorb.load()

