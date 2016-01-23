# ChemEvol
Python package to read in a star formation history file, input galaxy parameters and run a chemical evolution model to determine the evolution of gas, metals and dust in galaxies.

Running the script following the instructions below will produce:

1. a results data file with galaxy parameters in (filename.dat with name provided
  by user from inits file)
1. a pop-up plot for looking at results quickly

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040).

If you use this code, please do cite the above papers.  The license is provided with this package.

## Requirements

### Python packages
- numpy
- astropy
- logger
- matplotlib

### input files needed
The code reads in a star formation history from a file called filename.sfh.  This needs to be in the form: time (yr), SFR (Msolar/yr).    An example is provided with this code `Milkyway.sfh` based on the SFH for the Milky Way in Yin et al 2009 (A & A, 505, 497).

### input data needed
The code requires a dictionary of parameters to feed in, these are set in main.py and can be changed to suit following the comments.

## Running the code
The code can be run when in the directory by either `python main.py` or by using the following example (note: requires a SFH file called Milkyway.sfh and delayed.sfh provided
  with this package).  For further details on the definition of the parameters please see comments in `main.py`.

```python
import functions as f
import data as d
from evolve import ChemModel
'''
Initialise the parameters of the model
'''
initialise your galaxy parameters here and choice of models
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


'''
Initialise arrays
'''

snrate = []
all_results = []
galaxies = []

'''
Call the functions and run the chemical evolution model for however many galaxies
(list in inits) provided by the user
'''

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

 # write out table to file based on 'name' identifier provided by users
 name = item['name']
 t.write(str(name+'.dat'), format='ascii', delimiter=' ')

 # write each run through inits to produce final array galaxies for playing with
 galaxies.append(params)

# eg: to get any list of entries for any dictionary name for quick look
# [g['name'] for g in galaxies]
# [g['mgas'] for g in galaxies]
