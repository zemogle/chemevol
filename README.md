# ChemEvol
[![Build Status](https://travis-ci.org/zemogle/chemevol.svg?branch=master)](https://travis-ci.org/zemogle/chemevol)

Python package to read in a star formation history file, input galaxy parameters and run a chemical evolution model to determine the evolution of gas, metals and dust in galaxies.

Running the script following the instructions below will produce:

1. a results data file with galaxy parameters in (filename.dat with name provided
  by user from inits file - see instructions below)

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040), and De Vis et al 2017b (MNRAS, 471, 1743), with the latest developments described in De Vis et al 2021	.

The version of this code used for De Vis et al 2017b is in release [v_de_vis2017](https://github.com/zemogle/chemevol/releases/tag/v_de_vis2017).
The master version contains the version described in De Vis et al 2021.

#
If you use this code, please do cite the above papers.  The license is provided with this package.


## Installation

Download this repository and do the following:
```
python setup.py install
```

## Requirements

### Python packages
- numpy
- astropy
- matplotlib [Not a strict requirement]

### Input files needed
The code reads in a star formation history or star formation efficiency file from a file called filename.sfh or filename.sfe respectively.  This needs to be in the form: time (yr), SFR/SFE (Msolar/yr), inflow rates (Msolar/yr) and in the directory where you are running the code.   Examples are provided with this code in the examples directory, e.g. `Milkyway.sfh` based on the SFH for the Milky Way in Yin et al 2009 (A & A, 505, 497) or any of the other .sfh and .sfe files (described in De Vis et al 2017,2020).

### Input data needed
The code also requires a dictionary of initial parameters. This can be done by providing a .json or .csv file and using the package installed, or by running the example python script provided with this package.

There are example data files in the folder `<Download Dir>/chemevol/chemevol/examples/` which show the correct format of each type of file.  Feel free to copy these example data files into the directory where you wish to run the code and follow the instructions.  A detailed breakdown of the parameter files and inputs needed are found at the bottom of this readme.
The most basic example of using the chemevol package would be:

```python
from chemevol import ChemModel
ch = ChemModel(**item)
snrate = ch.supernova_rate()
all_results = ch.gas_metal_dust_mass(snrate)
```

There are helper functions for loading batch files in CSV and valid JSON format.
*Note*: Valid JSON uses double quotes for definitions and lower case for booleans, e.g. `"myvalue" : false`

## Running the code as a package using json or csv data files

```python
import chemevol as ch
galaxies = ch.BulkEvolve('<File directory>/data.json')
galaxies.upload_json()
galaxies.evolve_all()
```
Or for CSV files:
```python
import chemevol as ch
galaxies = ch.BulkEvolve('data.csv')
galaxies.upload_csv()
galaxies.evolve_all()
```
See the files in `<Download Dir>/chemevol/chemevol/examples/` for the correct format of each type of file (.py, .json, .csv examples are provided).

### Viewing the results
Once the code is run you will have an array called `galaxies.results` with all the parameters in.  To look at this data try:
```python
[g['dust_all'] for g in galaxies.results] #will print out the dust_all ( which is the total dust mass that was produced, without taking into account any dust destruction).
[g['mgas'] for g in galaxies.results] #will print out all the gas masses
gasmass  = galaxies.results[0]['mgas'] #etc
```

The code writes data to a file (if you use the example in `<Download Dir>/chemevol/chemevol/examples/data.json` the code writes out two files called `Model_A.dat` and `Model_B.dat`).  To read in this data for plotting or playing with you can use `astropy.table`:
```python
from astropy.table import Table
import matplotlib.pyplot as plt

t = Table.read('Model_A.dat', format='ascii')
plt.semilogy(t['fg'],t['dustmass']/(t['mgas']+t['mstars']))
```

Alternatively you can run the code without the bulkevolve class above:

#### What follows is a short description of the various scripts in the chemevol package and their general usage. ####

lookups.py, functions.py and evolve.py form the core of the chemevol package.
lookups.py contains the bulk of the hardcoded tables, such as the metal yield tables etc. 
When one wants to add different yield tables, these have to included in lookups.py in the appropriate format.

functions.py contains various functions that give numerical descriptions of the various formula used throughout the model calculations.
E.g. it contains all the dust processing mechanisms or the prescription for how the inflows/outflows and recycling is handled.

evolve.py runs the actual model, by keeping track of the gas, dust and metal content through each time-step as a result of all the ongoing physical processes.
To run a model, the appropriate inputs have to be defined and passed as a dictionary to the initialiser of the ChemModel class.

Next to these core files, there are some script in the examples directory that automate some of the setup and running of the chemical evolution models.
These call to the chemevol package and provide the easiest way to run one or more models.

example_multi.py gives some example input dictionaries, as well as the appropriate way to process them with the chemevol package and write the output data to .dat or .csv files.
Copy the example_multi.py file to the desired directory (where your .sfh file is and where you want to have your results). Then edit the parameters in the init dictionaries within this script.

run_grid_parallel.py reads input parameters from a predefined csv file and processes them in parallel. This script has all the functionallity to efficiently run and save a large amount of models. 

make_gridTable_dust.py creates the csv file for a grid of models, typically processed by run_grid_parallel.py.

make_sfe_sfh.py shows how to create SFE/SFH files with some examples. Any changes to the star formation history or the inflow rates will be made here.

There is also a log file that lists the warnings. These are warnings about corrections made during the model calculations. Actual errors or failures of the code will still cause the code to crash and give an appropriate error message.

**.sfe   star formation efficiency files, as defined in De Vis et al. (2020)
**.sfh   star formation history files, as defined in De Vis et al. (2017)
**.dat 	 some example models

## Full list of Input Parameters to Model

- name: 						someway to identify your galaxy
- gasmass_init: 					initial gas mass in solar masses
- starmass_init: 				initial stellar mass in solar masses
- dustmass_init: 				initial dust mass in solar masses
- Z_init: 						initial metallicity in solar masses
- SFH: 						filename for star formation history file (filename.sfh)
							if you don't specify a file, it will default to MW like SFH
- add_bursts:					if you are using a .sfe file, you can add bursts following De Vis et al 2020						
- t_end:						end of time array for the model
- t_start:						start of time array for the model
- gamma: 						power law for extrapolation of SFH if using SFH generated by MAGPHYS code, otherwise set to zero
- IMF_fn:          					choice of IMF function: Chab/chab/c, TopChab/topchab/tc, Kroup/kroup/k or Salp/salp/s
- dust_source: 					choice of dust sources to be included:
								SN: supernova dust only
								LIMS: low intermediate mass stars dust only
								LIMS+SN: both SN and LIMS included
								ALL: includes supernovae, LIMS and grain growth combined
- cold_gas_fraction:				Fraction of the gas that is in cold dense clouds, default is 0.5 eg Asano et al 2013
- use_THEMIS:					If set to True, use the THEMIS dust production and destruction mechanisms					
- delta_lims_fresh: 				Efficiency of fresh metals condensing into dust grains in LIMS (1<M_i<8Msun)
							Set to 0.15-0.4 in Morgan & Edmunds (2003); 0.15 in De Vis et al 2017b in press
- reduce_sn_dust				reduce the contribution from SN dust (currently set to values from
							Todini & Ferrera 2001).  Can be True or False. To reduce dust mass
							then quote number to reduce by (factor).
- destroy: 						on: add dust destruction from SN shocks: True or False
								mass: Amount of material destroyed by each SN
								(typically 1000 or 100Msun)
- fragmentgrains:				THEMIS photofragmentation rate of large a-C:H/a-C grains
- effective_snrate_factor: 			factor that the sn_rate needs to be corrected by to account for previous SN clearing dust in the vicinity of this SN, default is 0.36.
- graingrowth = 					grain growth parameter from THEMIS (clouds) or Mattsson & Andersen 2012 (depending on whether use_THEMIS is True)
- graingrowth2 = 				grain growth parameter from THEMIS (diffuse ISM)
- inflows: 						there are five parameters
								on: do you wish to turn inflows on: input True or False
								mass: 	If 0, keep the masses defined in sfh file. 
								If >0, set up the inflows so that `mass' gives the total amount of inflowing material
								inflows_metals = metallicity of inflow Y or N: input a number appropriate for primordial 
								inflow gas eg 1e-3 to 1e-4 (Rubin et al 2012, Peng & Maiolino 2013).
								inflows_xSFR = inflow rate of gas is X * SFR: input a number X; this parameter is not used in the current version, but can easily be 											reintroduced using the inflows_SFR function in functions.py
								inflows_dust = amount of dust inflow: input a number appropriate for dust eg 0.4 x the metallicity (Edmunds 2000)
- outflows: 					there are four parameters: Outflow rates are taken from Nelson et al (2019)
								on: do you wish to turn outflows on Y or N: input True or False
 								outflows_metals = metallicity of inflow: input True or False
						  		(True = metallicity of system, False = 0)
						 		outflows_dust = amount of dust in outflow: input True of False
								(True = dust/gas of system, False = 0)
								reduce = reduces the outflow rates from Nelson et al (2019) by this factor
- recycle: 						there are three parameters: 
								on: do you wish to turn outflow recycling on Y or N: input True or False
								esc_prob_perGyr: probability per Gyr that the outflowing material is lost to the IGM before it is recycled
								reaccr_time_factor: scale the recycling time up or down by this factor
- available_metal_fraction:			fraction of metals that is available for grain growth in the diffuse ISM, 
								this is effectively the maximum dust-to-metal ratio in the diffuse ISM.
								This is also used to calculate the maximum dust-to-metal ratio in clouds as 2.45 times the available_metal_fraction.
- SNyield:						String identifier to choose which SN metal yield table is used (see lookups.py).			
- AGByield:					String identifier to choose which SN metal yield table is used (see lookups.py).			
- totyields:						Boolean deciding whether total metal yields or fresh metal yields should ne used.				
- isotopes:						Isotopes/metal budgets to be tracked throughout the model. 
							These have to be consistent with the SNyields and AGByields tables (see lookups.py).			
- Pristine_isotope_fractions:		If there are metals in the inflows, Pristine_isotope_fractions decides how much of which metals are present.
