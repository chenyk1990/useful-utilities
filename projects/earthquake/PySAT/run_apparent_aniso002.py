import pysat
import matplotlib.pyplot as plt
import numpy as np
gammaArr = np.array([]); aziArr = np.array([]); phiArr = np.array([])

# eTensor=pysat.elasticTensor()
# eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
# eTensor.rot_dip_strike2(dip=33., strike=300.)
# eTensor.decompose_MN()

for strike in np.arange(0, 360, 1):
    print strike
    eTensor=pysat.elasticTensor()
    eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
    eTensor.rot_dip_strike2(dip=30., strike=strike)
    eTensor.decompose_MN()
    gamma, eps, delta = eTensor.etETI.get_thomsen()
    aziA, phifa =  eTensor.get_aziA()
    gammaArr=np.append(gammaArr, gamma)
    aziArr=np.append(aziArr, aziA)
    phiArr=np.append(phiArr, phifa)

eTensor=pysat.elasticTensor()
eTensor.set_radial(vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790)
gamma, eps, delta = eTensor.get_thomsen()
eTensor.rot_dip_strike2(dip=90., strike=80.)
aziA, phifa =  eTensor.get_aziA()

gammaArr= gammaArr/gammaArr.max()
aziArr= aziArr/aziArr.max()

strikeArr = np.arange(0, 360, 1)
plt.plot(strikeArr, gammaArr, 'r-', lw=5,label='radial anisotropy')
plt.plot(strikeArr, aziArr, 'b-', lw=5,label='azimuthal anisotropy')
plt.ylabel('normalized anisotropy', fontsize=30)
plt.xlabel('strike (deg)', fontsize=30)
plt.ylim([0.99, 1.01])
plt.legend()
plt.figure()

plt.plot(strikeArr[~ np.isnan(phiArr)], phiArr[~ np.isnan(phiArr)], 'b-', lw=5)
plt.ylabel('fast axis angle(deg)', fontsize=30)
plt.xlabel('strike (deg)', fontsize=30)
# plt.ylim([phiArr[~ np.isnan(phiArr)].mean()-0.1, phiArr[~ np.isnan(phiArr)].mean()+0.1])
plt.show()
# plt.plot(dipArr, gammaArr)




    

