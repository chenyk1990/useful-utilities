import pysat
import matplotlib.pyplot as plt
import numpy as np
gammaArr = np.array([]); aziArr = np.array([]); phiArr = np.array([])

eTensor=pysat.elasticTensor()
eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
eTensor.rot_dip_strike2(dip=33., strike=300.)
eTensor.decompose_MN()

for dip in np.arange(0, 90, 1):
    print dip
    eTensor=pysat.elasticTensor()
    eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
    eTensor.rot_dip_strike2(dip=dip, strike=30.)
    eTensor.decompose_MN()
    gamma, eps, delta = eTensor.etETI.get_thomsen()
    aziA, phifa =  eTensor.get_aziA()
    gammaArr=np.append(gammaArr, gamma)
    aziArr=np.append(aziArr, aziA)
    phiArr=np.append(phiArr, phifa)

eTensor=pysat.elasticTensor()
eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
gamma, eps, delta = eTensor.get_thomsen()
eTensor.rot_dip_strike2(dip=90., strike=30.)
aziA, phifa =  eTensor.get_aziA()

gammaArr= gammaArr
aziArr= aziArr

dipArr = np.arange(0, 90, 1)
plt.plot(dipArr, gammaArr/gammaArr.max(), 'r-', lw=5,label='radial anisotropy')
plt.plot(dipArr, aziArr/aziArr.max(), 'b-', lw=5,label='azimuthal anisotropy')
plt.ylabel('anisotropy', fontsize=30)
plt.xlabel('dip (deg)', fontsize=30)

plt.legend(loc='lower left')
# plt.figure()
# 
# plt.plot(dipArr[~ np.isnan(phiArr)], phiArr[~ np.isnan(phiArr)], 'b-', lw=5)
# plt.ylim([phiArr[~ np.isnan(phiArr)].mean()-0.1, phiArr[~ np.isnan(phiArr)].mean()+0.1])
# plt.ylabel('fast axis angle(deg)', fontsize=30)
# plt.xlabel('dip (deg)', fontsize=30)
plt.show()
# plt.plot(dipArr, gammaArr)




    

