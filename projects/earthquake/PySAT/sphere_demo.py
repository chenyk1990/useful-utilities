import pysat
import numpy as np
# 
eTensor=pysat.elasticTensor()
eTensor.elastic_DB(mtype='post-perovskite')

kceq=pysat.Christoffel(etensor=eTensor)
# kceq.sphere(dtheta=1., dphi=1., outfname='sphere005.asdf', group=True)
kceq.read_asdf(infname='sphere005.asdf')
# kceq.plot3d(ptype='rel', stype='rel')
kceq.plot2d(ptype='abs', stype='abs', ds=15, cmap='jet', fastpolar=False)
