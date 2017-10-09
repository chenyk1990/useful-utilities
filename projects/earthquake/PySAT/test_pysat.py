import pysat
import numpy as np
# 
eTensor=pysat.elasticTensor()
eTensor.elastic_DB(mtype='ice')
# eTensor.set_love(A=15.55, C=11.96, L=4.76, N=5.33, F=3.99, mtype='HTI')
# eTensor.Voigt2Cijkl()
eTensor.rot_dip_strike(dip=30., strike=50.)

eTensor2=pysat.elasticTensor()
eTensor2.elastic_DB(mtype='ice')
# eTensor.set_love(A=15.55, C=11.96, L=4.76, N=5.33, F=3.99, mtype='HTI')
# eTensor.Voigt2Cijkl()
eTensor2.rot_dip_strike(dip=30., strike=50., method='euler')