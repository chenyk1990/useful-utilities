import pysat
import numpy as np

for dip in np.arange(0, 90, 1):
    for strike in np.arange(0, 360, 1):
        
        eTensor=pysat.elasticTensor()
        eTensor.set_radial(vsv=3.57, vsh=4.0, vpv=6.14, vph=7.0, eta=0.7, rho=2790)
        # eTensor.elastic_DB(mtype='ol')
        eTensor.rot_dip_strike2(dip=dip, strike=strike)
        eTensor2=pysat.elasticTensor()
        eTensor2.set_radial(vsv=3.57, vsh=4.0, vpv=6.14, vph=7.0, eta=0.7, rho=2790)
        # eTensor2.elastic_DB(mtype='ol')
        g = pysat.euler2mat(dip, strike+270., 0, 'syzx')
        eTensor2.Cijkl=np.einsum('ia,jb,kc,ld,abcd->ijkl', g, g, g, g, eTensor2.Cijkl)
        eTensor2.Cijkl2Voigt()
        
        print np.allclose(eTensor.Cijkl, eTensor2.Cijkl)
        

