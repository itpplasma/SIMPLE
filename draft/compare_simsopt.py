#%%
import simsopt
import simsopt.mhd
import booz_xform

wout_file = "wout.nc"

b = booz_xform.Booz_xform()
b.read_wout(wout_file)
b.run()
# %%

v = simsopt.mhd.vmec.Vmec(wout_file)

# %%
