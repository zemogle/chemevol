'''
This script generates a csv file for a grid of models (which is generally processed by run_grid_parallel.py).
The grid values are set in the steps list for each parameter defined in params.
Additional parameters can be added if special care is taken to insert them in the correct place in the string to be written to the csv.
When making changes, it is key to ensure the indices of mastertable are still correct. 
'''


import numpy as np
import math
from subprocess import call
import random


'''
To make a grid of models, we here specify which parameters are varied, and what values should be used within the grid.
The example parameter values are from the grid of dust models in De Vis et al (2020)
'''
params=["initial_mass","SFH",'include_bursts',"IMF",'fc','LIMS_yield',"reduced sn", "destruction", "fragmentation", "grain growth cloud", "grain growth diffuse" ,'reaccr_time_factor','available_metals','SNyield','AGByield']
#steps=[["1e7","3e7","1e8","3e8","1e9","3e9","1e10","3e10","1e11","3e11","1e12"],["average.sfe","fast.sfe","slow.sfe"],[False],["Chab","Salp"],[0.5],[0.15],[5],[15],[0.5],[4000],[5],[2.,4.],[0.3],["tot_LC18_R150"],['KA18_high']]
steps=[["1e7","3e7","1e8","3e8","1e9","3e9","1e10","3e10","1e11","3e11"],["average.sfe","fast.sfe","slow.sfe"],[True,False],["Chab","Salp","TopChab","Kroup"],[0.5],[0.15],[5],[15],[0.05],[4000],[5],[2.,1.,0.5,0.25,0.1],[0.2],["tot_LC18_R000","tot_LC18_R150","tot_LC18_R300","MA92_ori_extra","MM02_000"],['Nugrid','FRUITY','KA18_low','KA18_high','KA10','VG97']]
# Next we initialise and then populate a table which contains all permutations off the above parametervalues.
mastertable=[]

def addParam(table,steps):
	if table==[]:
		newtable=steps
	else:	
		newtable=list(np.zeros(len(table)*len(steps)))
		for ii in range(len(table)):
			for i in range(len(steps)):
				newtableline=np.append(table[ii],steps[i])
				newtable[ii*len(steps)+i]=newtableline
	return newtable	

def generateName(table,i):
	name="GridZb_%s"%(i+1)
	return name #[0:-1]


for i in range(len(params)):
	mastertable=addParam(mastertable,steps[i])



'''
Once the mastertable with all permutations is created, we generate the csv file by including 
each of the values in each line of the mastertable into the full input string in the csv.
The other values are set to a default value.
If additional parameters are to be varried, they need to be inserted in the string at the right place (following the header of the csv).
'''
outtable=open("modelGridZ.csv","w")   		
outtable.write("# name, gasmass_init, starmass_init, dustmass_init, Z_init, SFH, add_bursts, t_end, t_start,\
 gamma, IMF_fn, dust_source, cold_gas_fraction, use_THEMIS, delta_lims_fresh, reduce_sn_dust_on, reduce_sn_dust_factor,\
 destroy_on, mass_destroy, fragment_on, fragment_tau, effective_snrate_factor, epsilon_grain, epsilon_grain2,\
 inflows_on, inflows_mass, inflows_metals, inflows_xSFR, inflows_dust, outflows_on,outflows_metals, outflows_dust, outflows_reduce,\
 recycle_on, esc_prob_perGyr, reaccr_time_factor, available_metal_fraction, SNyield, AGByield,totyields,isotopes,Pristine_isotope_fractions \n")
for i in range(len(mastertable)):
	dustsource='ALL' 
	if np.float(mastertable[i][10])==0 and np.float(mastertable[i][9])==0:
		dustsource='SN+LIMS'
	name=generateName(mastertable[i],i)	
	totyields=False	
	if "tot_" in mastertable[i][13]:
		totyields=True
		mastertable[i][13]=mastertable[i][13][4::]	
	dust_reduce=True
	if mastertable[i][6]==0:
		dust_reduce=False	
	dust_destroy=True
	if float(mastertable[i][7])==0:
		dust_destroy=False
	dust_frag=True
	if float(mastertable[i][8])==0:
		dust_frag=False	
	string="%s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s, %s,  %s,  %s, %s,  %s,  %s,  %s,  %s,  %s,  %s, %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s, %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s,  %s \n"%(name,\
		mastertable[i][0], 0., 0., 0., mastertable[i][1], mastertable[i][2], 13.8, 1., 0, mastertable[i][3], dustsource, mastertable[i][4], True, mastertable[i][5],  dust_reduce, mastertable[i][6],\
		dust_destroy, mastertable[i][7], dust_frag, mastertable[i][8],0.36, mastertable[i][9], mastertable[i][10], \
		True, 0., 0., 1., 0., True, True, True, 1., True, 0.2, mastertable[i][11], mastertable[i][12], mastertable[i][13], mastertable[i][14], totyields, 'Z;O;N', '1.;0.435;0.055')
	print(string)
	outtable.write(string)
outtable.close()
	
