# ChemEvol
[![Build Status](https://travis-ci.org/zemogle/chemevol.svg?branch=master)](https://travis-ci.org/zemogle/chemevol)

Python package to read in a star formation history file, input galaxy parameters and run a chemical evolution model to determine the evolution of gas, metals and dust in galaxies.

Running the script following the instructions below will produce:

1. a results data file with galaxy parameters in (filename.dat with name provided
  by user from inits file - see instructions below)

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040), with latest features discussed in De Vis et al (submitted 2016).

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
The code reads in a star formation history from a file called filename.sfh.  This needs to be in the form: time (yr), SFR (Msolar/yr) and in the directory where you are running the code.   An example is provided with this code `Milkyway.sfh` based on the SFH for the Milky Way in Yin et al 2009 (A & A, 505, 497).

### Input data needed
The code also requires a dictionary of initial parameters. These can be set directly by adding dictionaries directly into main.py, or providing a CSV or JSON file.  There are example files in the folder `<Download Dir/chemevol/chemevol/examples/` which show the correct format of each type of file.  Feel free to copy these example data files into the directory where you wish to run the code.

There are helper functions for loading batch files in CSV and valid JSON format.
*Note*: Valid JSON uses double quotes for definitions and lower case for booleans, e.g. `"myvalue" : false`

## Running the code

```python
import chemevol as ch
galaxies = ch.BulkEvolve('<File directory>/data.json')
galaxies.upload_json()
galaxies.evolve_all()
```
Or for CSV files:
```python
import chemevol as ch
galaxies = ch.BulkEvolve('<File directory>/data.csv')
galaxies.upload_csv()
galaxies.evolve_all()
```
See the files in `<Download Dir>/chemevol/chemevol/examples/` for the correct format of each type of file.

### Viewing the results
Once the code is run you will have an array called `galaxies.results` with all the parameters in.  To look at this data try:
```python
[g['dust_all'] for g in galaxies.results] #will print out the dust_all
[g['mgas'] for g in galaxies.results] #will print out all the gasmasses
gasmass  = galaxies.results[0]['mgas'] #etc
```

The code writes data to a file (if you use the example in `<Download Dir>/chemevol/chemevol/examples/data.json` the code writes out two files called `Model_I.dat` and `Model_VI.dat`).  To read in this data for plotting or playing with you can use `astropy.table`:
```python
from astropy.table import Table
import matplotlib.pyplot as plt

t = Table.read('Model_VI.dat', format='ascii')
plt.semilogy(t['fg'],t['dustmass']/(t['mgas']+t['mstars']))
```
