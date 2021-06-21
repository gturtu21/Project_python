import main_input_MOD 

path='/home/giorgio/dpd_metrangolo/gaff_charges/fdg3/FD3_IN_WATER_ETHANOL_5PC/'
trajfile = 'din.nc'
parmfile = 'fdg3_wt.prmtop'

mynewtraj=my_traj(path+trajfile,path+parmfile, 'ETN', 2)

mynewtraj.all_analysis()



