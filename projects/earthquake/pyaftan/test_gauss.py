import obspy
import pyaftan
import matplotlib.pyplot as plt
import numpy as np
# 
# prefname='/lustre/janus_scratch/life9360/PRE_PHP/ses3d_2016_R/E000.98S43.pre'
# f='SES.174S110.SAC'
# f='./sac_data/SES.98S43.SAC'
# f='/work3/leon/COR_TEST/2008.APR/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
# f='/work3/leon/COR_TEST/2009.NOV/COR/109C/COR_109C_BHZ_R21A_BHZ.SAC'
f='/work3/leon/COR_TEST/COR/TA.109C/COR_TA.109C_BHZ_TA.R21A_BHZ.SAC'
st=obspy.read(f)
tr=st[0]
tr1=pyaftan.aftantrace(tr.data, tr.stats)
tr1.makesym()
tr1.write('o.sac', format='sac')
fdata=tr1.gaussian_filter_aftan(0.05)
tr1.data=fdata
tr1.write('filtered.sac', format='sac')