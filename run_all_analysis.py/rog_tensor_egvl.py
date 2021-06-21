import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import pytraj as pt
file_tensor=pt.read_pickle('GYR_TENSOR_JAN1_COS_TFN_NDEN_2.pk')

### file_tensor is a tuple, the second elemtn of the tuple is the tensor itself

tensor=file_tensor[1]

len_traj=len(tensor)

eigenvalues=np.zeros((3,len_traj))

for i in range(len_traj):
    F=np.array([[tensor[i,0], tensor[i,3], tensor[i,4]],
            [tensor[i,3], tensor[i,1], tensor[i,5]],
            [tensor[i,4], tensor[i,5], tensor[i,2]]])
    eigenvalues[:,i]=LA.eigh(F)[0]

zeta1=np.transpose(eigenvalues[0,:])
zeta2=np.transpose(eigenvalues[1,:])
zeta3=np.transpose(eigenvalues[2,:])
componenti=np.column_stack((zeta1, zeta2, zeta3))
print(np.shape(componenti))
plt.plot(zeta1,label='zeta1')
plt.plot(zeta2,label='zeta2')
plt.plot(zeta3,label='zeta3')
plt.xlabel('Number of Steps')
plt.ylabel(R'Principal Moments of Gyration Tensor($\AA^2)$')
plt.legend()
plt.show()
np.savetxt('tensor_components.txt',componenti,fmt='%.6f')
