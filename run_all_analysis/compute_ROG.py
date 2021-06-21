import pytraj as pt
import numpy as np
from numpy import linalg as LA


def radius_of_gyration(traj,compute_tensor=True):
    """ compute the tensor of gyration of
    the JAN molecules """
    for solute in range(traj.ndendr):
        if compute_tensor:
            #pt.to_pickle(pt.radgyr(traj, ':' + str(solute+1)+' tensor',dtype='dict'),'GYR_TENSOR_JAN'+str(solute+1)+'.pk')
            # pt.radgyr_tensor returns a tuple of len 2:
            # - the first element of the tuple is the radius of gyration over time:
            # - the second element is a tuple of 6-dimension arrays, these 6 elements are the components
            #   XX , YY , ZZ , XY, XZ, YZ of the gyration tensor
            #pt.to_pickle(pt.radgyr_tensor(traj, ':' + str(solute+1)),'GYR_TENSOR_JAN'+str(solute+1)+'_COS_'+traj.cosolvent+'_NDEN_'+str(traj.ndendr)+'.pk')
            results = pt.radgyr_tensor(traj, ':' + str(solute+1))
            tensor = results[1]
            len_traj = len(tensor)
            eigenvalues = np.zeros((3,len_traj))

            for i in range(len_traj):
                F=np.array([[tensor[i,0], tensor[i,3], tensor[i,4]],
                        [tensor[i,3], tensor[i,1], tensor[i,5]],
                        [tensor[i,4], tensor[i,5], tensor[i,2]]])
                eigenvalues[:,i]=LA.eigh(F)[0]
                
            eigenvalues=np.transpose(eigenvalues)
            np.savetxt('tensor_components_JAN'+str(solute+1)+'_COS_'+traj.cosolvent+'_NDEN_'+str(traj.ndendr)+'.txt',eigenvalues,fmt='%.6f')


        else:
            pt.to_pickle(pt.radgyr(traj, ':' + str(solute+1)),'ROG_JAN'+str(solute+1)+'.pk')
    return 




