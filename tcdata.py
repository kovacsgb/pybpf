"""
Module containing a useable interface for the data output of Budapest-Florida code
in Python.

# Current features

- History class providing support for fort.18 files
- Model class providing support for fort.95 files
- LimmitCycle class providing support for fort.19 files
- bpfDataRead function provides easy initialization in given directory

Created by: Gábor B. Kovács 2021
"""
from astropy.constants import cgs
import numpy as np
import astropy.constants as constants
from astropy import units


class BaseData:
    """
    Base class for the three most common Bupest-Florida Code output, such as
     - fort.95 which contains the static starting model profile
     - fort.19 which contains the profiles for two full periods
     - fort.18 which contains the photosphere data, for the entire run
    """
    def __init__(self, data=None, col_names=None):
        self.datablock=data
        self.column_names=col_names

    def data(self,key):
        if key in self.datablock:
            return self.datablock[key]
        else:
            raise KeyError('There is no'+str(key)+'data member.')
    
    def __getattr__(self,key):
        return self.data(key)
    

class History(BaseData):
    """
    This class provides standard interface for working with fort.18 files.
    Column names:
     - time:    Elapsed time from hydro start in days
     - vel_rad: Radial velocity in the photosphere
     - L_phot:  Photosphere luminosity
     - radius:  Radius of surface (One zone below nzn.)
     - L_surf:  Luminosity of the surface zone
    """
    columnnames=('time','vel_rad','L_phot','radius','L_surf') 
    def __init__(self,path):
        data,colnames=self.read_data(path)
        super().__init__(data,colnames)
        

    def read_data(self,path):
            datafile=open(path)
            dataLines=datafile.readlines()
#This should be static!
            columns=[]
            for i in range(0,len(dataLines[0].split())):
                columns.append([float(x.split()[i]) for x in dataLines])
                columns[i]=np.array(columns[i])

            return dict(zip(self.columnnames,columns)), self.columnnames

            
class Model(BaseData):
    """
    Class providing interface for working with fort.95 files.
    """
    columnnames=('zone','radius','pressure','spec_vol','energy','temperature','turbulent_energy',
                    'dm','L_c','p_turb','L_r','sk','L_turb','s','cs','fcsl','st','cla','taup','taud','opacity',
                    'edt','cp')
    def __init__(self,path):
        data,colnames=self.read_data(path)
        super().__init__(data,colnames)
        
    def read_data(self,path):
        datafile=open(path)
        dataLines=datafile.readlines()
        columns=[]
        for i in range(0,len(dataLines[0].split())-1):
            columns.append([float(x.split()[i].replace('D','E')) for x in dataLines])
            columns[i]=np.array(columns[i])
        return dict(zip(self.columnnames,columns)), self.columnnames


class DataPoint(BaseData):
    """
    Helper class for fort.19 data access, holds one record of fort.19 data.
    Fields are:
    phase       - phase of pulsation
    zone        - zone number of mass cell
    velocity    - velocity of mass cell
    L_r         - radiative luminosity
    radius      - radius at mass cell
    temperature - temperature of the cell
    dm          - size of the cell
    opacity     - opacity in the cell
    energy      - internal energy
    e_t         - turbulent energy
    entropy     - specific entropy
    F_t         - turbulent luminosity
    F_c         - convective luminosity
    pressure    - gas pressure in cell
    p_t         - turbulent pressure
    p_eddy      - viscous pressure of the cell
    mass        - mass (Langrangian) coordinate of the cell.
    """
    columnnames=('phase','zone','velocity','L_r','radius','temperature','dm','opacity','energy','e_t','entropy',
                'F_t','F_c','pressure','p_t','p_eddy','mass')

    def __init__(self,Line : str):
        #TODO: d D e E csere!

        data=dict(zip(DataPoint.columnnames,[float(x.translate({"d": "e","D" : "e"})) for x in Line.split()]))
        #for x in Line.split(): print(float(x))
        super().__init__(data,DataPoint.columnnames)

    @classmethod
    def createDataPointList(cls):
        return dict(zip(DataPoint.columnnames, [ list() for i in enumerate(cls.columnnames)]))
    
    @classmethod
    def fillDataPointList(cls,datapoints : list):
        if not all(isinstance(x,cls) for x in datapoints):
            raise TypeError("datapoints should be a list of DataPoint")
        datapoint_list=cls.createDataPointList()
        for point in datapoints:
            for key in cls.columnnames:
                datapoint_list[key].append(point.data(key))
        for key in datapoint_list:
                datapoint_list[key]=np.array(datapoint_list[key])
        return datapoint_list


class RawProfiles:
    """
    Container holding DataPoints extracted from fort.19
    """
    def __init__(self, filename):
        self.datablock=None
        self.num_time_series=0
        self.num_profiles=0
        self.read_file(filename)

    def read_file(self,filename):
        infile=open(filename)
        lines = infile.readlines()
        self.datablock = [DataPoint(line) for line in lines]
        self.num_time_series = int(self.datablock[-1].zone)
        self.num_profiles = int(self.datablock[-1].phase)
        infile.close()

    def __getitem__(self,key) -> DataPoint:
        return self.datablock[key]



class TimeSeries(BaseData):
    """
    Holding one zones time series of the limit cycle
    """
    def __init__(self,rawprofile_obj : RawProfiles,zone_number):
        if not all([isinstance(rawprofile_obj,RawProfiles),isinstance(zone_number,(int,float))]):
            raise TypeError("Wrong arguments")
        super().__init__(self.createTimeSeries(rawprofile_obj,zone_number),DataPoint.columnnames)

    def createTimeSeries(self,rawprofile_obj : RawProfiles, zone_number) -> dict:
        datapoints=rawprofile_obj[zone_number::rawprofile_obj.num_time_series]
        #print(type(datapoints),zone_number,rawprofile_obj.num_time_series,len(rawprofile_obj.datablock))
        data=DataPoint.fillDataPointList(datapoints)
        return data

class Profiles(BaseData):
    """
    Holding a limit cycle profile at given phase number.
    """
    def __init__(self,rawprofile_obj : RawProfiles,phase_number):
        super().__init__(self.createProfile(rawprofile_obj,phase_number),DataPoint.columnnames)
    
    def createProfile(self,rawprofile_obj : RawProfiles,phase_number):
        datapoints=[point for point in rawprofile_obj if point.phase == phase_number]
        data=DataPoint.fillDataPointList(datapoints)
        return data

class LimitCycle:
    """
    Class providing user interface for the fort.19 data, both timeseries and profiles extraction, using the same interface.
    For columnnames see DataPoint class.
    """
    columnames=DataPoint.columnnames
    def __init__(self,infile_name):
        if isinstance(infile_name,str):
            rawdata=RawProfiles(infile_name)
        elif isinstance(infile_name,RawProfiles):
            rawdata = infile_name
        else:
            raise TypeError("Parameter should be either string or RawProfile object.")
        self.timeSeries=[TimeSeries(rawdata,zone_number) for zone_number in range(1,rawdata.num_time_series+1)]
        self.profiles=[Profiles(rawdata,phase_number) for phase_number in range(1,rawdata.num_profiles+1)]
        self.num_time_series = len(self.timeSeries)
        self.num_profiles=len(self.profiles)


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
    if "return_rawprofile" in kwargs:
        return_rawprofile = kwargs["return_rawprofile"]

    path_to_fort95 = path+"/fort.95"
    path_to_fort18 = path+"/fort.18"
    path_to_fort19 = path+"/fort.19"

    model_obj=Model(path_to_fort95)
    history_obj=History(path_to_fort18)
    rawprofile_obj=RawProfiles(path_to_fort19)
    limitcycle_obj=LimitCycle(rawprofile_obj)
    if return_rawprofile:
        return model_obj,history_obj,limitcycle_obj,rawprofile_obj
    else:
        return model_obj,history_obj,limitcycle_obj




#testing utilities
if(__name__ == '__main__'):
    from matplotlib import pyplot as plt
    import math
    a_point = DataPoint("45\t136\t12\t50\t5.6\t2e5\t4.13e+29\t125.68\t1.25e40\t1.2e38\t12000\t1e29\t12e5\t1.5e4\t1e2\t0.98\n")
    fort19_path='/home/gabesz/SPHERLS_playground/bpf_konv/konv-b/fort.19'#'/home/gabesz/Dokumentumok/VS_Code_workspace/calcion/fort.19'
    fort19_data=RawProfiles(fort19_path)
    print(a_point.datablock)
    print(a_point.dm)
    print(fort19_data[15].zone)

    #Testing Timeseries:
    #surface:
    surface_element=fort19_data.num_time_series-2
    print("érték:",fort19_data.num_time_series)
    surface_series=TimeSeries(fort19_data,surface_element)
    #print(surface_series.radius)
    mylog=np.vectorize(math.log10)
    #plt.plot(surface_series.phase/77.,surface_series.F_r/constants.L_sun.cgs)
    
    aProfile=Profiles(fort19_data,1)

    fort19_handler=LimitCycle(fort19_path)

    for phase in range(0):#,len(fort19_handler.profiles)):
        #phase=40
        ts_xdata=fort19_handler.timeSeries[fort19_handler.num_time_series-2].phase
        ts_ydata=fort19_handler.timeSeries[fort19_handler.num_time_series-2].L_r
        xdata=fort19_handler.profiles[phase].zone
        ydata=fort19_handler.profiles[phase].F_c
        #print(xdata,ydata)
        plt.subplot(2,1,1)
        plt.plot(xdata,ydata/constants.L_sun.cgs)
        plt.plot(xdata,fort19_handler.profiles[phase].L_r/constants.L_sun.cgs)
        plt.plot(xdata,fort19_handler.profiles[phase].F_t/constants.L_sun.cgs)
        #plt.plot(aProfile.zone[1:146],mylog(aProfile.temperature[1:146]))
        plt.subplot(2,1,2)
        plt.plot(ts_xdata/77,ts_ydata/constants.L_sun.cgs)
        plt.scatter(ts_xdata[phase]/77,ts_ydata[phase]/constants.L_sun.cgs,color="red")
        plt.show()
        plt.pause(5)


    #Testing reader
    mod,his,lim = bpfDataRead("/home/gabesz/SPHERLS_playground/bpf_konv/konv-b/")

    #Other tests
    print(len(lim.profiles), lim.num_profiles)
