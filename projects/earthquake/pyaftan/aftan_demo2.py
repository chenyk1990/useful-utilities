import pyaftan
import obspy
import matplotlib.pyplot as plt
import numpy as np



tr=obspy.read('./test.sac')[0]
atr1=pyaftan.aftantrace(tr.data, tr.stats)
atr2=pyaftan.aftantrace(tr.data, tr.stats)
# aftan analysis using pyaftan
atr1.aftan(tmin=2., tmax=40., phvelname='ak135.disp')
# aftan analysis using compiled fortran library
atr2.aftanf77(tmin=2., tmax=40., phvelname='ak135.disp')
# compare two aftan results (no phase matched filter)
atr1.ftanparam.FTANcomp(atr2.ftanparam, compflag=2 )
plt.show()
# compare two aftan results (with phase matched filter)
atr1.ftanparam.FTANcomp(atr2.ftanparam, compflag=4 )
plt.show()