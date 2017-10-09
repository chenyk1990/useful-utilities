import pysat
import numpy as np
# 
eTensor=pysat.elasticTensor()
# eTensor.elastic_DB(mtype='ol')
eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)

kceq=pysat.Christoffel(etensor=eTensor)
# kceq.set_direction_cartesian(pv=[1.1,2.2,3.3])
kceq.set_direction_cartesian(pv=[1., 0., 0.])
kceq.get_phvel()
# kceq.get_grad_mat()
# kceq.get_group_velocity()