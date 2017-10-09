import pysat
import numpy as np
# 
eTensor=pysat.elasticTensor()
eTensor.elastic_DB(mtype='ol')
# eTensor.set_love(A=15.55, C=11.96, L=4.76, N=5.33, F=3.99, mtype='HTI')
# 
# eTensor.Cvoigt =  np.array([[ 47.601,   44.371, 46.342, 0,    0., 0.],
#                             [  44.371,  18.906, 20.327,  0,    0., 0.],
#                             [ 46.342,  20.327, 18.839,    0.0,    0., 0.],
#                             [0,  0,   0.0,  3.1112,    0., 0.],
#                             [    0.,     0.,   0.0,     0., 5.218, 0],
#                             [    0.,     0.,   0.0,     0.,0,  4.5369]])
# eTensor.Voigt2Cijkl()
eTensor.rotT([0,0,1], 30.)
 
# eTensor.rot_dip_strike(dip=43., strike=120., method='default')
eTensor2=pysat.elasticTensor()
eTensor2.elastic_DB(mtype='ol')
# eTensor2.set_love(A=15.55, C=11.96, L=4.76, N=5.33, F=3.99, mtype='HTI')
# eTensor2.Cvoigt[:] =  np.array([[ 47.601,   44.371, 46.342, 0,    0., 0.],
#                                         [  44.371,  18.906, 20.327,  0,    0., 0.],
#                                         [ 46.342,  20.327, 18.839,    0.0,    0., 0.],
#                                         [0,  0,   0.0,  3.1112,    0., 0.],
#                                         [    0.,     0.,   0.0,     0., 5.218, 0],
#                                         [    0.,     0.,   0.0,     0.,0,  4.5369]])
eTensor2.rotB([0,0,1], 30.)

# eTensor3=pysat.elasticTensor()
# eTensor3.elastic_DB(mtype='ice')

# eTensor2.rot_dip_strike(dip=43., strike=120., method='euler')


# R=pysat.euler2mat(10,20,30)
# eTensor.Voigt2Cijkl()
# eTensor.set_radial(8.7, 4.5, 3.3, 0.1, 0.5, 1)
# eTensor.set_love(A=13.961, C=15.013, L=3.21, N=3.404, F=5.765)
# eTensor.set_thomsen(vp=5.3, vs=3.2, eps=0.1, gamma=0.05, delta=0.00, rho=4000)
# # eTensor.Voigt2Cijkl()
# eTensor.rotB([0,1,0], 90.)
# eTensor.rotB([1,0,0], 90.)
# eTensor.rotTB([0,0,1], 90.)
# eTensor.Cijkl2Voigt()
# 
# # 
# eTensor2=pysat.elasticTensor()
# # # eTensor.Voigt2Cijkl()
# # 
# eTensor2.set_love(A=13.961, C=15.013, L=3.21, N=3.404, F=5.765)
# # eTensor.Voigt2Cijkl()
# eTensor2.rotT([0,1,0], 90.)

# eTensor2.rotB([1,0,0], 90.)
# eTensor2.rotB([0,0,1], 90.)