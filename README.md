# chemevol
Chemical evolution python package

# Requirements

## Python packages
- numpy
- astropy
- logger
- matplotlib

## input files needed
The code reads in a star formation history from a file called filename.sfh.  This needs to be in the form time (yr), SFR (Msolar/yr).    An example is provided `MilkyWay.sfh` based on the SFH for the Milky Way in Yin et al 2009 (A & A, 505, 497).

## input data needed
The code requires a dictionary of parameters to feed in, these are set in main.py and can be changed to suit.

## Running the code
The code can be run using the following example along with a SFH file.  

```python
import plots as p
from evolve import ChemModel
inits = {
        'gasmass_init': 4.8e10,
				'SFH': 'MilkyWay.sfh',
        't_end': 20.,
				'gamma': 0,
				'IMF_fn': 'Chab',
				'dust_source': 'ALL',
				'destroy': True,
				'inflows':{'metals': 0., 'xSFR': 0, 'dust': True},
				'outflows':{'metals': True, 'xSFR': 0, 'dust': True},
				'cold_gas_fraction': 0.5,
				'epsilon_grain': 1e5,
        'destruct': 1000.
              }
ch = ChemModel(**inits)
time, mgas = ch.gas_mass()
time, mstars = ch.stellar_mass()
time, metalmass, metallicity = ch.metal_mass(mgas)
snrate = ch.supernova_rate()
time, mdust, dust_sources, dust_metals, timescales = ch.dust_mass(mgas,metallicity,snrate)

gasfraction = mgas/(mgas+mstars)

p.figure(time,mgas,mstars,metalmass,metallicity,mdust,dust_metals,gasfraction,dust_sources,timescales)
```
