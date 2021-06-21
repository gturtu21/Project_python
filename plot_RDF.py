import pickle 
import matplotlib.pyplot as plt
import numpy as np

a=int(input('traj 1 or 5?'))

if a == 1:
    with open('novolume_rdf_traj1.pickle','rb') as pickle_file:
    	content = pickle.load(pickle_file)

elif a ==5: 
    with open('novolume_rdf_traj5.pickle','rb') as pickle_file:
        content = pickle.load(pickle_file)


couples=list(content.keys())
print(list(content.keys()))
print(content[couples[0]][1])

for el in couples:
    np.savetxt('rdf_'+el+'.txt', content[el][1],fmt='%.6f')


############################################################################
#fix, axs = plt.subplots(2)
#axs[0].plot(content[':JAN@%f:TFN@%ho'][0],content[':JAN@%f:TFN@%ho'][1])
#axs[1].plot(content2[':JAN@%f:TFN@%ho'][0],content2[':JAN@%f:TFN@%ho'][1])
#plt.show()

