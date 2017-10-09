import obspy
import pyaftan
import matplotlib.pyplot as plt
import numpy as np
# 
# prefname='/lustre/janus_scratch/life9360/PRE_PHP/ses3d_2016_R/E000.98S43.pre'
# f='SES.174S110.SAC'
# f='./sac_data/SES.98S43.SAC'
f='/work3/leon/COR_TEST/2008.APR/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
# f='/work3/leon/COR_TEST/2009.NOV/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
# f='/work3/leon/COR_TEST/COR/TA.109C/COR_TA.109C_BHZ_TA.R21A_BHZ.SAC'
st=obspy.read(f)
tr=st[0]
tr1=pyaftan.aftantrace(tr.data, tr.stats)
# tr1.stats.sac.b=-0.0
tr1.makesym()
tr1.aftanf77(piover4=-1., pmf=True, vmin=2.4, vmax=3.5, ffact=1. , tmin=4.0, tmax=30.0)
tr1.get_snr()
d1=tr1.ftanparam.arr2_2[8,:tr1.ftanparam.nfout2_2]

# f='/work3/leon/COR_TEST/2008.APR/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
f='/work3/leon/COR_TEST/2009.NOV/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
# f='/work3/leon/COR_TEST/COR/TA.109C/COR_TA.109C_BHZ_TA.R21A_BHZ.SAC'
st=obspy.read(f)
tr=st[0]
tr1=pyaftan.aftantrace(tr.data, tr.stats)
# tr1.stats.sac.b=-0.0
tr1.makesym()
tr1.aftanf77(piover4=-1., pmf=True, vmin=2.4, vmax=3.5, ffact=1. , tmin=4.0, tmax=30.0)
tr1.get_snr()
d2=tr1.ftanparam.arr2_2[8,:tr1.ftanparam.nfout2_2]

# f='/work3/leon/COR_TEST/2008.APR/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
# f='/work3/leon/COR_TEST/2009.NOV/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
f='/work3/leon/COR_TEST/COR/TA.109C/COR_TA.109C_BHZ_TA.R21A_BHZ.SAC'
st=obspy.read(f)
tr=st[0]
tr1=pyaftan.aftantrace(tr.data, tr.stats)
# tr1.stats.sac.b=-0.0
tr1.makesym()
tr1.aftanf77(piover4=-1., pmf=True, vmin=2.4, vmax=3.5, ffact=1. , tmin=4.0, tmax=30.0)
tr1.get_snr()
d3=tr1.ftanparam.arr2_2[8,:tr1.ftanparam.nfout2_2]
