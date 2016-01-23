'''
Chemevol - Python package to read in a star formation history file,
input galaxy parameters and run a chemical evolution model to determine the evolution
of gas, metals and dust in galaxies.

Running this script will produce
(a) a results data file (filename.dat) with file name given by user
(b) a pop-up plot for looking at gas, dust and metal evolution

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040).

If you use this code, please do cite the above papers.

Copyright (C) 2015 Haley Gomez, Edward Gomez and Simon Schofield, Cardiff University and LCOGT
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

********************************************************************************
'''

'''------------------------------------------------------------------------
First set up initial parameters for galaxy model by editing the dictionary
initial_galaxy_params

- name: 				someway to identify your galaxy
- gasmass_init: 		initial gas mass in solar masses
- SFH: 					filename for star formation history file (filename.sfh)
						if you don't specify a file, it will default to MW like SFH
- t_end:				end of time array for chemical integrals
- gamma: 				power law for extrapolation of SFH
						if using SFH generated by MAGPHYS code
- IMF_fn:          		choice of IMF function: Chab/chab/c, TopChab/topchab/tc,
		  				Kroup/kroup/k or Salp/salp/s
- dust_source: 			choice of dust sources to be included:
						SN: supernova dust only
						LIMS: low intermediate mass stars dust only
						LIMS+SN: both SN and LIMS included
						ALL: includes supernovae, LIMS and grain growth combined
- reduce_sn_dust		reduce the contribution from SN dust (currently set to values from
						Todini & Ferrera 2001).  If leave default specify False. To reduce dust mass
						then quote number to reduce by
- destroy: 				add dust destruction from SN shocks: True or False
- inflows: 				there are two parameters
 						metals = metallicity of inflow: input a number
								 xSFR = inflow rate is X * SFR: input a number
								 dust = amount of dust inflow: input a number
- outflows: 			there are two parameters
 						metals = metallicity of inflow: input True or False
						True = metallicity of system, False = 0
								 xSFR = inflow rate is X * SFR: input a number
								 dust = amount of dust in outflow: input True of False
							  	 		True = dust/gas of system, False = 0
- cold_gas_fraction = 	fraction of gas in cold dense state for grain growth
					  	typically 0.5-0.9 for high z systems, default is 0.5
- epsilon_grain = 		grain growth parameter from Mattsson & Andersen 2012
						default is 500 for t_grow ~ 10Myr.
- destruct = 			amount of material destroyed by each SN
						(typically 1000 or 100Msun)


Each run will be used to generate the evolution of dust, gas,
SFR, metals and stars over time
---------------------------------------------------------------------------
'''

import functions as f
from evolve import ChemModel
import data as d
from astropy.table import Table
import matplotlib.pyplot as plt

# initialise your galaxy parameters here and choice of models
# each {} entry is per galaxy separated by comma in list

inits = [
			{	'name': 'Model_I',
				'gasmass_init': 4e10,
				'SFH': 'Milkyway.sfh',
        		't_end': 20.,
				'gamma': 0,
				'IMF_fn': 'Chab',
				'dust_source': 'ALL',
				'reduce_sn_dust': False,
				'destroy': False,
				'inflows':{'metals': 0., 'xSFR': 0, 'dust': 0},
				'outflows':{'metals': False, 'xSFR': 0, 'dust': False},
				'cold_gas_fraction': 0.5,
				'epsilon_grain': 800,
        		'destruct': 0  },

			{	'name' : 'Model_IV',
				'gasmass_init': 4e10,
  				'SFH': 'delayed.sfh',
          		't_end': 20.,
  				'gamma': 0,
  				'IMF_fn': 'Chab',
  				'dust_source': 'ALL',
  				'reduce_sn_dust': 6,
  				'destroy': True,
  				'inflows':{'metals': 0., 'xSFR': 1.5, 'dust': 0},
  				'outflows':{'metals': True, 'xSFR': 1.5, 'dust': True},
  				'cold_gas_fraction': 0.5,
  				'epsilon_grain': 800,
          		'destruct': 100.  }
        ]

snrate = []
all_results = []
galaxies = []

for item in inits:
	ch = ChemModel(**item)
	'''
	call modules to run the model:
	snrate: 		SN rate at each time step - this also sets time array
					so ch.supernova_rate() must be called first to set
					time array for the entire code

	all results: 	t, mg, m*, mz, Z, md, md/mz, sfr,
					dust_source(all), dust_source(stars),
					dust_source(ism), destruction_time, graingrowth_time
	'''
	snrate = ch.supernova_rate()
	all_results = ch.gas_metal_dust_mass(snrate)
	# write all the parameters to a dictionary for each init set
	params = {'time' : all_results[:,0],
		   'mgas' : all_results[:,1],
		   'mstars' : all_results[:,2],
		   'metalmass' : all_results[:,3],
		   'metallicity' : all_results[:,4],
		   'dustmass' : all_results[:,5],
		   'dust_metals_ratio' : all_results[:,6],
		   'sfr' : all_results[:,7],
		   'dust_all' : all_results[:,8],
		   'dust_stars' : all_results[:,9],
		   'dust_ism' : all_results[:,10],
		   'time_destroy' : all_results[:,11],
		   'time_gg' : all_results[:,12]}
	params['fg'] = params['mgas']/(params['mgas']+params['mstars'])
	params['ssfr'] = params['sfr']/params['mgas']
	# write to astropy table
	t = Table(params)
	# write out to file based on 'name' identifier
	name = item['name']
	t.write(str(name+'.dat'), format='ascii', delimiter=' ')
	# if you want an array including every inits entry:
	galaxies.append(params)

# make some quick look up plots
#d.figure(time,mgas,mstars,metalmass,metallicity,dustmass,dust_metals_ratio,gasfraction)
