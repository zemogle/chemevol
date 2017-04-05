'''
Example of using ChemEvol without using the BulkEvolve class, but with multiple galaxies.

---------------------------------------------------
First set up initial parameters for galaxy model by editing the dictionary
initial_galaxy_params

- name: 				someway to identify your galaxy
- gasmass_init: 		initial gas mass in solar masses
- SFH: 					filename for star formation history file (filename.sfh)
						if you don't specify a file, it will default to MW like SFH
- t_end:				end of time array for chemical integrals
- gamma: 				power law for extrapolation of SFH if using SFH generated by MAGPHYS code, otherwise set to zero
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
- destroy: 				on: add dust destruction from SN shocks: True or False
						mass: Amount of material destroyed by each SN
						(typically 1000 or 100Msun)
- inflows: 				there are four parameters
							on: do you wish to turn inflows on: input True or False
							inflows_metals = metallicity of inflow Y or N: input a number appropriate for primordial
							inflow gas eg 1e-3 to 1e-4 (Rubin et al 2012, Peng & Maiolino 2013).
							inflows_xSFR = inflow rate of gas is X * SFR: input a number X
							inflows_dust = amount of dust inflow: input a number appropriate
							for dust eg 0.4 x the metallicity (Edmunds 2000)
- outflows: 			there are two parameters
							on: do you wish to turn outflows on Y or N: input True or False
 							outflows_metals = metallicity of inflow: input True or False
						   	(True = metallicity of system, False = 0)
						 	outflows_dust = amount of dust in outflow: input True of False
							(True = dust/gas of system, False = 0)
							(The equation for outflows is take from simulations of feedback in Hopkins et al 2012, see Eq 27 in Feldmann et al 2015)
- cold_gas_fraction = 	fraction of gas in cold dense state for grain growth
					  	typically 0.5-0.9 for high z systems, default is 0.5 eg Asano et al 2013
- epsilon_grain = 		grain growth parameter from Mattsson & Andersen 2012
						default is 500 for t_grow ~ 10Myr.

Each run will be used to generate the evolution of dust, gas,
SFR, metals and stars over time
---------------------------------------------------------------------------
'''

from chemevol import ChemModel
from astropy.table import Table

# initialise your galaxy parameters here and choice of models
# each {} entry is per galaxy separated by comma in list

inits = [
			{	'name': 'Model_A',
				'gasmass_init': 4e10,
				'SFH': 'Milkyway.sfh',
        		't_end': 20.,
				'gamma': 0,
				'IMF_fn': 'Chab',
				'dust_source': 'All',
				'reduce_sn_dust': False,
				'destroy': {'on': False, 'mass': 0},
				'inflows':{'on': False, 'metals': 0., 'xSFR': 0, 'dust': 0},
				'outflows':{'on': False, 'metals': False, 'dust': False},
				'cold_gas_fraction': 0.5,
				'epsilon_grain': 0},

			{	'name': 'Model_B',
				'gasmass_init': 4e10,
				'SFH': 'MilkyWay.sfh',
        		't_end': 20.,
				'gamma': 0,
				'IMF_fn': 'Chab',
				'dust_source': 'All',
				'reduce_sn_dust': False,
				'destroy': {'on': True, 'mass': 0},
				'inflows':{'on': False, 'metals': 0., 'xSFR': 0, 'dust': 0},
				'outflows':{'on': False, 'metals': False, 'dust': False},
				'cold_gas_fraction': 0.5,
				'epsilon_grain': 800}
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
