from .tcdata import *
from .calcion import *

def bpfDataRead(path : str, **kwargs):
    """
    Function to read in a Budapest-Florida-Code working directory, return a tuple of objects of Model, History and LimitCycle
    classes. The only parameter is the relative/absolute path to the BpF working directory.
    Additionally it returns a RawProfile object too.
    Other keywords:
    return_rawprofile -> if true it returns the RawProfile object too
    """
    #TODO: It should be refactorized if it will be extended for other datafiles too.
    return_rawprofile = False
    do_ionization = False
    X = Y = 0.
    if "return_rawprofile" in kwargs:
        return_rawprofile = kwargs["return_rawprofile"]
    if "do_ionization" in kwargs:
        if 'X' in kwargs and 'Y' in kwargs:
            do_ionization = kwargs["do_ionization"]
            X=kwargs["X"]
            Y=kwargs['Y']
        else:
            raise RuntimeError("No metallicity was given")
        
    path_to_fort95 = path+"/fort.95"
    path_to_fort18 = path+"/fort.18"
    path_to_fort19 = path+"/fort.19"

    model_obj=Model(path_to_fort95)
    history_obj=History(path_to_fort18)
    rawprofile_obj=RawProfiles(path_to_fort19)
    if (do_ionization):
        rawprofile_obj = Ionization.IonizationForRawprofile(rawprofile_obj,X,Y)
    limitcycle_obj=LimitCycle(rawprofile_obj)
    if return_rawprofile:
        return model_obj,history_obj,limitcycle_obj,rawprofile_obj
    else:
        return model_obj,history_obj,limitcycle_obj


