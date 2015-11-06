# chemevol
Chemical evolution python package

# Requirements

## Python packages
- numpy
- astropy

## input files needed
The code reads in a star formation history from a file called filename.sfh.  This needs to be in the form time (yr), SFR (Msolar/yr).    An example is provided `MilkyWay_yin.sfh` based on Yin et al 2009 (A & A, 505, 497).

## Running the code
```python
from evolve import ChemModel
inits = {
          'gasmass_init':4e10,
          'SFH':'MilkyWay.sfh',
          'gamma':0,
          'IMF_fn':'Chab',
          'dust_source':'ALL',
          'destroy':True,
          'inflows':0,
          'outflows':0
          }
ch = ChemModel(**inits)
ch.gas_mass(choice='chab')
```
