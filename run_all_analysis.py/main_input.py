import os, sys
import numpy as np
import pytraj as pt
import matplotlib.pyplot as plt
import pickle 

############
print(sys.path)
try:
    sys.path.remove('/usr/local/amber16/lib/python2.7/site-packages')
except ImportError:
    print('mpi4py has the wrong path!!')
    pass
print(sys.path)

#####################



################################# IMPORT THE TRAJECTORY ##################################
home=os.getenv('HOME')
print(home)

which_path=int(input('Which MD trajectory to analyze? \n \
        - 1: Two fdg3 in etoh 5%; \n \
        - 2: One fdg3 in etoh 5% ; \n \
        - 3: One fdg3 in water; \n \
        - 4: Two fdg3 in fetoh (0-100 ns) ; \n \
        - 5: Two fdg3 in fetoh (100-200 ns) ; \n \
        - 6: Two fdg3 in fetoh (0-200 ns); \n \
        - 7: Two fdg3 in etoh (0-200 ns); \n \
Print the number correspoding to the desired trajectory = '))
if which_path == 1:
    path ='/dpd_metrangolo/gaff_charges/fdg3/box_water_ethanol_5pc_100/production_merged/'
    trajfile = 'dintot.nc'
    parmfile = 'dimer_wtet5pc.prmtop'
if which_path == 2:
    path='/dpd_metrangolo/gaff_charges/fdg3/FD3_IN_WATER_ETHANOL_5PC/'
    trajfile = 'din.nc'
    parmfile = 'fdg3_wtet5pc.prmtop'
if which_path == 3:
    path = '/dpd_metrangolo/gaff_charges/fdg3/box_just_water/production/'
    trajfile = 'din.nc'
    parmfile = 'fdg3_wt.prmtop'
if which_path == 4:
    path = '/dpd_metrangolo/gaff_charges/fdg3/box_water_tfethanol_5pc_100_close_dendrons/production/' 
    trajfile = 'din.nc'
    parmfile = 'dimer_wtfet5pc.prmtop'
if which_path == 5:   ### insert ethanol_tainah
    path = '/dpd_metrangolo/gaff_charges/fdg3/box_water_tfethanol_5pc_100_close_dendrons/production/' 
    trajfile = 'din_rst.nc' 
    parmfile = 'dimer_wtfet5pc.prmtop'
if which_path == 6:
    path = '/dpd_metrangolo/gaff_charges/fdg3/box_water_tfethanol_5pc_100_close_dendrons/production/'
    trajfile = 'dintot.nc'
    parmfile = 'dimer_wtfet5pc.prmtop'
if which_path == 7:
    path ='/dpd_metrangolo/gaff_charges/fdg3/box_water_ethanol_5pc_100/production_tainah/'
    trajfile= 'din.nc'
    parmfile= 'dimer_wtet5pc.prmtop'

full_path=home+path
all_traj = pt.iterload(full_path+trajfile,full_path+parmfile,[(0,-2)])
##########################################################################################
###### check files belongs to the proper class ######

class my_traj(pt.TrajectoryIterator):
    """ cosolvent should be 'ETN' or 'TFN' or None """
    def __init__(self,filename, top, cosolvent, number_of_dendrons, timerange=[(0,-2)]):
        self.cosolvent = cosolvent
        self.ndendr=number_of_dendrons
        self.time = timerange
        super().__init__(filename,top)
        pt.iterload(filename,top,self.time)
    def radius_of_gyration(self):

        #pt.radgyr(self,':JAN')
        pass
    def fluorine_atomnames(self):
        return [atom.name for atom in self.top.atoms if atom.type=='f' if atom.resid==0]
    def compute_rdf(self, save_to_file=True, plot_results=True):
        from compute_RDF import all_rdf 
        results = all_rdf(self) 
        if save_to_file:
            import pickle
            with open('novolume_rdf_'+ mynewtraj.cosolvent+'.pickle', 'wb') as handle:
                pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        if plot_results:
            import matplotlib.pyplot as plt
            for couples in results.keys():
                plt.plot(results[couples][0], results[couples][1])
                plt.show()
        return self


#### load the trajectory as my_traj object ########################################
if which_path in (1,7):
    mynewtraj=my_traj(full_path+trajfile,full_path+parmfile, 'ETN', 2)
elif which_path == 2:
    mynewtraj=my_traj(full_path+trajfile,full_path+parmfile, 'ETN', 1)
    mynewtraj=mynewtraj['!:WAT']
elif which_path == 3:
    mynewtraj=my_traj(full_path+trajfile,full_path+parmfile, None , 1)
elif which_path in (4,5,6) :
    mynewtraj=my_traj(full_path+trajfile,full_path+parmfile, 'TFN' , 2)


####################################################################################

#print(mynewtraj.fluorine_atomnames())

import time
start=time.time()
###############  THE FOLLOWING ANALYSIS WILL BE CARRIED OUT: ######

which_analysis=input('Which analysis: (choose RDF or DIHEDRAL or DISTCH or DISTFF or ROG or ALL)')


# 1) RDF;

if which_analysis == 'RDF':
    from compute_RDF import all_rdf 
    results = all_rdf(mynewtraj)
    ## SAVE DICTIONARY TO FILE 
    with open('novolume_rdf_traj'+str(which_path)+'.pickle', 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

##########################################################################

# 2) EXTRACT DIHEDRAL

elif which_analysis == 'DIHEDRAL':
    from compute_dihedral import triazol_dihedral
    import matplotlib.pyplot as plt
    results = triazol_dihedral(mynewtraj)
    np.savetxt('dihedrals_media_traj'+str(which_path)+'.txt',results[1],fmt='%.3f')
    dihedrals = np.transpose(results[0])
    np.savetxt('dihedrals_traj'+str(which_path)+'.txt', dihedrals ,fmt='%.3f')

# 3) EXTRACT 'C' - 'ho' distances

elif which_analysis == 'DISTCH':
    from compute_distances import CH_distances
    results = CH_distances(mynewtraj)
    for key in results[0]:
        plt.plot(results[1][key],label=key)
    plt.legend()
    plt.show()


# 4) EXTRACT 'f' - 'f' intermolecular distances

elif which_analysis == 'DISTFF':
    start = time.time()
    from extract_distance import interm_distances
    results = interm_distances(mynewtraj)
    end = time.time()
    print('Total time to extract fluorine intermolecular distances : ' + str(end-start) + 'min')

elif which_analysis == 'ROG':
    start = time.time()
    from compute_ROG import radius_of_gyration
    radius_of_gyration(mynewtraj,True)

end=time.time()
print(end-start)
