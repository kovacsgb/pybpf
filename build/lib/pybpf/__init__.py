"""
PyBpf
=====

This package is meant for a modern interface for dealing with the
outputs of Budapest-Florida code (Yecko et al. 1996.).

Current Features
----------------

 - Full interface for fort.95,fort.18 and fort.19 files
 - OOP datastructure for easy acces all the data
 - Ionization fraction calculations for the fort.19 data
 - A reader function to easily initialize all data from a BpF code-
 directory

Usage
-----

You can import it as:

>>> import pybpf

and you can easily get data from full bpf working directory (assumed
that fort.95 fort.19 and fort.18 files DOES exists.)

>>> working_path = 'path/to/working/dir'
>>> mod, hist, lim = pybpf.bpfDataRead(working_path,\
     do__ionization = True, X=0.75,Y=0.2496)

and you have initialized the whole directory.
"""


from . import tcdata
from . import calcion
from .bpfreader import bpfDataRead
