import transformations
import math
import transforms3d
R1 = transformations.rotation_matrix(math.pi/2, [0, 1, 0])
R2 = transforms3d.axangles.axangle2mat([0, 1, 0], math.pi/2)
R2 = transforms3d.axangles.axangle2mat([-1/3., 2/3., 2/3.], math.pi/180.*-74)