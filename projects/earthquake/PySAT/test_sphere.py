import pysat
import numpy as np
# 
eTensor=pysat.elasticTensor()
# eTensor.elastic_DB(mtype='ice')
eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
eTensor.rot_dip_strike2(dip=20., strike=300.)
eTensor.info='Hexagonal Symmetric Media (lower hemisphere)'


kceq=pysat.Christoffel(etensor=eTensor)
# kceq.sphere(dtheta=1., dphi=1., outfname='jiayi002.asdf', group=True)
# # kceq.sphere(dtheta=5., dphi=5., outfname='sphere002.asdf', group=True)
kceq.read_asdf(infname='jiayi002.asdf')
# kceq.plot3d(ptype='abs', stype='abs')
kceq.plot2d(ptype='abs', stype='abs', ds=10, theta0=180., hsph='lower', cmap='cv')
# kceq.set_direction_cartesian(pv=[1.1,2.2,3.3])
# kceq.set_direction_cartesian(pv=[1.1,2.2,3.3])
# kceq.get_phvel()
# kceq.get_grad_mat()
# kceq.get_group_velocity()