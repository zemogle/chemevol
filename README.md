# ChemEvol
Python package to read in a star formation history file, input galaxy parameters and run a chemical evolution model to determine the evolution of gas, metals and dust in galaxies.

Running the script following the instructions below will produce:

1. a results data file with galaxy parameters
1. a pop-up plot for looking at results quickly

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040).

If you use this code, please do cite the above papers.

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
The code can be run when in the directory using the following example (note: requires a SFH file).  

```python
import functions as f
import data as d
from evolve import ChemModel

inits = {
        'gasmass_init': 4e10,
        'SFH': 'Milkyway.sfh',
        't_end': 20.,
        'gamma': 0,
        'IMF_fn': 'Chab',
        'dust_source': 'ALL',
        'destroy': True,
        'inflows':{'metals': 0., 'xSFR': 0, 'dust': 0},
        'outflows':{'metals': True, 'xSFR': 1.5, 'dust': True},
        'cold_gas_fraction': 0.5,
        'epsilon_grain': 1000.,
        'destruct': 1000.
        }

ch = ChemModel(**inits)

snrate = ch.supernova_rate()

dust_sources, timescales, all_results = ch.gas_metal_dust_mass(snrate)

time = all_results[:,0]
mgas = all_results[:,1]
mstars = all_results[:,2]
metalmass = all_results[:,3]
metallicity = all_results[:,4]
mdust = all_results[:,5]
dust_metals = all_results[:,6]
sfr = all_results[:,7]

gasfraction = mgas/(mgas+mstars)
ssfr = sfr/mstars

d.writedata(time, mgas, mstars, sfr, ssfr, mdust, metalmass, metallicity, gasfraction)

d.figure(time,mgas,mstars,metalmass,metallicity,mdust,dust_metals,gasfraction,dust_sources,timescales)

```
