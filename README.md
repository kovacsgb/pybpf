PyBpf package
=============

This package is meant for a modern interface for dealing with the
outputs of Budapest-Florida code (Yecko et al. 1996.).

Current Features
----------------

 - Full interface for fort.95,fort.18 and fort.19 files
 - OOP datastructure for easy access all the data
 - Ionization fraction calculations for the fort.19 files
 - A reader function to easily initialize all data from a BpF code-directory

Prerequisities
------------

- python 3.4 or higher
- numpy
- astropy

Installation
------------

Use the setup.py script as superuser:

```bash
python3 setup.py install
```

Usage
-----

You can import it as:

```python
>>> import pybpf
```
and you can easily get data from full bpf working directory (assumed
that fort.95 fort.19 and fort.18 files DOES exists.)

```python
>>> working_path = 'path/to/working/dir'
>>> mod, hist, lim = pybpf.bpfDataRead(working_path,\
     do__ionization = True, X=0.75,Y=0.2496)
```
and you have initialized the whole directory.

Main data structure is in `pybpf.tcdata`:

 1. `pybpf.tcdata.Model` holds all data from fort.95 files, which is the
    static modell output of `tcmodel`
 2. `pybpf.tcdata.History` holds data from fort.18 files, created during
    `tcrelax` calculation descibing the state of the photosphered during
    the whole hydro run.
 3. `pybpf.tcdata.RawProfiles` holds data of fort.19 files in list format.
 4. `pybpf.tcdata.LimitCycle` holds data of fort.19 files, which describe
    the whole model in `isn/ipst` timepoint during the last two periods.
    It is extended with specific volume data. And can be extend by `calcion`
    module

`pybpf.calcion` Module can add ionization calculation for `RawProfiles`
object using Newton-iteration to solve the Saha-equations.

With the above function call we have the three object from the working directory,
calculated ionization fractions for H and He. We can any time print out the column
names of a given object:

```python
>>> print(mod.column_names)
('zone', 'radius', 'pressure', 'spec_vol', 'energy', 'temperature', 'turbulent_energy', 'dm', 'L_c', 'p_turb', 'L_r', 'sk', 'L_turb', 's', 'cs', 'fcsl', 'st', 'cla', 'taup', 'taud', 'opacity', 'edt', 'cp')
>>> print(hist.column_names)
('time', 'vel_rad', 'L_phot', 'radius', 'L_surf')
>>> print(lim.profile[0].column_names)
('HII_fraction', 'spec_vol', 'phase', 'energy', 'pressure', 'dm', 'opacity', 'mass', 'radius', 'zone', 'entropy', 'HeII_fraction', 'e_t', 'HeIII_fraction', 'p_t', 'L_r', 'F_c', 'F_t', 'temperature', 'p_eddy', 'velocity')
```
`LimitCycle` objects contains the time series for every zone and the profiles for every timepoint. We can use the fields `timeSeries` and `profiles` to access these (in a list format).
To access a datacolumn use `.` and the column name.

For example the light curve of the star can be accessed:

```python
>>> import matplotlib.pyplot as plt
>>> plt.plot(lim.timeSeries[-2].time/lim.num_profiles*2,lim.timeSeries[-2].L_r)
>>> plt.show()
```
Licence
-------
Distributed under MIT licence, for details see  LICENCE
