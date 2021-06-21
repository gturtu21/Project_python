import pytraj as pt
import numpy as np



def CH_distances(traj, time_averaged=True, save_to_file=True, plot=True):
    """ return a python array of all the distances between the Carbon-sp3 holding 
    the fluorinated structure (AtomName: 'C') and the 8 hydroxyl groups represented by
    alchoholic hydrogen (AtomType: 'ho') as a function of simulation timesteps """
    import numpy as np

    if traj.ndendr > 1:
        my_hydrogen_names = [atoms.name for atoms in traj.top.atoms if atoms.type=='ho']
        for hydrogen_i in my_hydrogen_names:
            certain_CH_couple = ':1@C ' + hydrogen_i
            pt.distance(traj, certain_CH_couple)
            #pt.pmap(pt.distance(traj, certain_CH_couple),traj)
    else:
        my_hydrogen_names = [atoms.name for atoms in traj.top.atoms if atoms.type=='ho']
        print(my_hydrogen_names)
        distances = {}
        for hydrogen_i in my_hydrogen_names:
            certain_CH_couple = ':1@C :1@' + hydrogen_i
            distances[hydrogen_i] = pt.distance(traj, certain_CH_couple)
    return my_hydrogen_names, distances


##########################################################################################################

def FF_distances(traj):
    """ return a python array of all the INTER-molecular distances between
    fluorine atoms in the dendrons present in the simulation box"""
    if traj.ndendr < 2:
        print('This method is available only for trajectory with more than 1 solute molecule')
        return 
    fluorine_dendron1=([atom for atom in traj.top.atoms if atom.type=='f' if atom.resid==0])
    fluorine_dendron2=([atom for atom in traj.top.atoms if atom.type=='f' if atom.resid==1])
    mymask = [':'+str(fluorines1.resid+1)+'@'+str(fluorines1.name)+' :'+str(fluorines2.resid+1)+'@'+str(fluorines2.name) \
    for fluorines1 in fluorine_dendron1 for fluorines2 in fluorine_dendron2]
    distances=pt.distance(traj, mymask)
    return distances
