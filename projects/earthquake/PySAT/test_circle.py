import pysat
import numpy as np
# 
eTensor1=pysat.elasticTensor()
# eTensor.elastic_DB(mtype='ice')
# eTensor1.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
eTensor1.set_radial(vsv=3.48, vsh=3.63, vpv=5.94, vph=6.28, eta=0.82, rho=2730)
gamma1, eps1, delta1 =  eTensor1.get_thomsen()
# eTensor1.set_thomsen(vp=6.5, vs=3.6, eps=0.1, gamma=0.1, delta=0.2, rho=2790)
# eTensor1.set_thomsen(vp=5.78, vs=3.46, eps=0.2, gamma=0.2, delta=0.21, rho=2720)
eTensor1.rot_dip_strike2(dip=34., strike=0.)
# eTensor.info='Hexagonal Symmetric Media (lower hemisphere)'


kceq1=pysat.Christoffel(etensor=eTensor1)
kceq1.circle(theta=90., dphi=1., group=True)
kceq1.plot_circle(showfig=False)


eTensor2=pysat.elasticTensor()
# eTensor2.set_radial(vsv=3.54, vsh=3.71, vpv=6.15, vph=6.47, eta=0.74, rho=2790)
eTensor2.set_radial(vsv=3.45, vsh=3.61, vpv=6.06, vph=6.24, eta=0.72, rho=2730)
gamma2, eps2, delta2 =  eTensor2.get_thomsen()
# eTensor2.set_thomsen(vp=6.5, vs=3.487, eps=0.1, gamma=0.1, delta=0.00, rho=2790)
# eTensor2.set_thomsen(vp=5.78, vs=3.414, eps=0.2, gamma=0.2, delta=-0.1, rho=2720)
eTensor2.rot_dip_strike2(dip=27., strike=90.)
# eTensor.info='Hexagonal Symmetric Media (lower hemisphere)'

kceq2=pysat.Christoffel(etensor=eTensor2)
kceq2.circle(theta=90., dphi=1.,  group=True)
kceq2.plot_circle(rel=True,polar=False, showfig=True)