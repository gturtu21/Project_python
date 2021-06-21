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

    def fluorine_atomnames(self):
        return [atom.name for atom in self.top.atoms if atom.type=='f' if atom.resid==0]

    #### AVAILABLE ANALYSIS IMPLEMENTED FOR THIS CLASS 

    def radius_of_gyration(self):
        from compute_RDF import all_rdf
        results = all_rdf(mynewtraj)
        with open('novolume_rdf_traj'+str(which_path)+'.pickle', 'wb') as handle:
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return 
    def extract_dihedral(self):
        from compute_DIHEDRAL import triazol_dihedral
        results = triazol_dihedral(self)
        np.savetxt('dihedrals_media_traj'+str(which_path)+'.txt',results[1],fmt='%.3f')
        dihedrals = np.transpose(results[0])
        np.savetxt('dihedrals_traj'+str(which_path)+'.txt', dihedrals ,fmt='%.3f')

    def extract_rdf(self, couples='all'):
        from compute_RDF import all_rdf
        if couples=='all':
            results=all_rdf(self)
            with open('novolume_rdf_traj'+str(which_path)+'.pickle', 'wb') as handle:
                pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return 

    def extract_ff_distances(self):
        from compute_DISTANCES import FF_distances
        distances = FF_distances(self)
        return distances

    def extract_c_ho_distances(self):
        from compute_DISTANCES import CH_distances
        distances = CH_distance(self)
        return distances

    ####### THIS METHOD PERFORMS ALL THE ANALYSIS LISTED ABOVE

    def all_analysis(self):
        import time
        start = time.time()
        self.radius_of_gyration()
        print('ROG analysis: Done')
        self.extract_dihedral()
        print('Dihedral analysis: Done')
        self.extract_rdf()
        print('RDF analysis: Done')
        self.extract_ff_distances()
        print('FF distances: Done')
        self.extract_c_ho_distances()
        print('C-ho distances: Done')
        end = time.time()
        print('Total time of execution: ', end-start//60 ,'minutes.')
        return 

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


