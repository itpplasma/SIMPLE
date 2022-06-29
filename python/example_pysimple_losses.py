#%%
from pysimple import simple, simple_main, params as p, new_vmec_stuff_mod as stuff
import numpy as np
import matplotlib.pyplot as plt

for item in dir(p):
    if item.startswith("__"): continue
    try:
        print(f'{item}: {getattr(p, item)}')
    except:
        print(f'{item}: NULL')

stuff.multharm = 3     # Fast but inaccurate splines
p.ntestpart = 32
p.trace_time = 1e-3
p.contr_pp = -1e10     # Trace all passing passing
p.startmode = -1       # Manual start conditions

tracy = p.Tracer()

simple.init_field(tracy, 'wout.nc',
    stuff.ns_s, stuff.ns_tp, stuff.multharm, p.integmode)

print(stuff.ns_s, stuff.ns_tp, stuff.multharm, p.integmode)

p.params_init()
#%%
# s, th_vmec, ph_vmec, v/v0, v_par/v
p.zstart = np.array([
  0.50000000000000000,        426.12215225438388,        1006.7214965304997,        1.0000000000000000,       0.54400086402893066,
  0.50000000000000000,        576.87315661554896,        1362.6901665761575,        1.0000000000000000,      -0.50067949295043945,
  0.50000000000000000,        69.461434003554956,        163.71928567067661,        1.0000000000000000,      -0.60121059417724609,
  0.50000000000000000,        119.76566838547805,        282.64080406435204,        1.0000000000000000,      -0.33900833129882812,
  0.50000000000000000,        1003.6136621639662,        2370.2605684794466,        1.0000000000000000,      -0.94397830963134766,
  0.50000000000000000,        1101.1103826203237,        2599.0543153399426,        1.0000000000000000,       0.32778775691986084,
  0.50000000000000000,        520.37070986619835,        1228.7216569466175,        1.0000000000000000,      -0.73738372325897217,
  0.50000000000000000,        537.03378741729728,        1267.7322734175166,        1.0000000000000000,       0.47406613826751709,
  0.50000000000000000,        1033.6784199821761,        2440.3712654285337,        1.0000000000000000,       0.61794042587280273,
  0.50000000000000000,        1315.7047615842671,        3105.9092888436548,        1.0000000000000000,       0.62871646881103516,
  0.50000000000000000,        421.72188303943193,        995.61271338657923,        1.0000000000000000,      -0.87669932842254639,
  0.50000000000000000,        28.736370014760556,        68.096062329317121,        1.0000000000000000,       0.43280315399169922,
  0.50000000000000000,        86.784471098573832,        205.11049049273603,        1.0000000000000000,       -2.8518438339233398E-002,
  0.50000000000000000,        1193.7735825402397,        2818.4268902570448,        1.0000000000000000,        9.7157359123229980E-002,
  0.50000000000000000,        377.18765979538455,        890.14823990452021,        1.0000000000000000,       0.52513253688812256,
  0.50000000000000000,        536.88778817668174,        1267.4043031984925,        1.0000000000000000,       0.26797866821289062,
  0.50000000000000000,        871.11512552510862,        2057.1528515801028,        1.0000000000000000,      -0.38690841197967529,
  0.50000000000000000,        1068.2010154596321,        2522.2447795671537,        1.0000000000000000,       0.11991834640502930,
  0.50000000000000000,        426.94545897383790,        1008.2662224117026,        1.0000000000000000,       0.92701447010040283,
  0.50000000000000000,        994.31400119800821,        2347.1102601822263,        1.0000000000000000,      -0.67407238483428955,
  0.50000000000000000,        407.05635936158160,        961.36752357409250,        1.0000000000000000,      -0.23453140258789062,
  0.50000000000000000,        336.24820278214833,        793.76381920314543,        1.0000000000000000,       0.63027143478393555,
  0.50000000000000000,        1258.5010114196741,        2970.7459030425866,        1.0000000000000000,       0.71755731105804443,
  0.50000000000000000,        1320.1347040907317,        3116.5802125732089,        1.0000000000000000,      -0.92375659942626953,
  0.50000000000000000,        39.516762021799494,        92.629824542685725,        1.0000000000000000,      -0.68713080883026123,
  0.50000000000000000,        711.72643833162078,        1680.0866099478583,        1.0000000000000000,       0.65814459323883057,
  0.50000000000000000,        241.50855519149621,        569.94103639470381,        1.0000000000000000,       0.10471045970916748,
  0.50000000000000000,        250.88311115728993,        592.77965868834599,        1.0000000000000000,       0.40150451660156250,
  0.50000000000000000,        471.91039660314226,        1113.8429071562314,        1.0000000000000000,      -0.49668657779693604,
  0.50000000000000000,        1184.1932075847819,        2795.8830683930264,        1.0000000000000000,        6.7270159721374512E-002,
  0.50000000000000000,        588.45543644634654,        1390.0327445620323,        1.0000000000000000,       0.40752923488616943,
  0.50000000000000000,        17.546104402881152,        41.910623029791218,        1.0000000000000000,       0.85118532180786133
]).reshape(32, 5).T
#%%
simple_main.run(tracy)

print(p.times_lost)

t = np.linspace(p.dtau/p.v0, p.trace_time, p.ntimstep)

plt.figure()
plt.semilogx(t, p.confpart_pass + p.confpart_trap)
plt.xlim([1e-4, p.trace_time])
plt.xlabel('time')
plt.ylabel('confined fraction')

plt.figure()
condi = np.logical_and(p.times_lost > 0, p.times_lost < p.trace_time)
plt.semilogx(p.times_lost[condi], p.perp_inv[condi], 'x')
plt.xlim([1e-4, p.trace_time])
plt.xlabel('loss time')
plt.ylabel('perpendicular invariant')
plt.show()

# %%