"""
Module containing all necessary calculations, functions and classes, for ionization calculations from
Budapest-Florida code output.

It uses the interface given in tcdata.py module.

Most important feature is:

Ionization.IonizationForRawprofile static method calculates ionization for given rawProfile object.

"""
#TODO:
# - Callable functions from outside
# - Testing utilites
# - resolve circular dependencies

from math import pi
from math import exp
from math import sqrt
import astropy.constants as consts
import astropy.units as u
from astropy.units.core import _recreate_irreducible_unit
import numpy as np
#from . import tcdata as tcdata
#import __init__
#import tcdata.tcdata as tcdata
try:
    from . import tcdata as tcdata
except (SystemError,ImportError):
    #this should be removed when all tests are defined in tests
    import tcdata as tcdata

class IterationError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class ConstantsCalculator:
    """
    Class calculating the right hand side of Saha equations.
    """
    R_HII = 0.0
    R_HeII = 0.0
    R_HeIII = 0.0
    n_H = 0.0
    n_He = 0.0
    M_He = 6.6464764e-24
    M_H = 1.6735575e-24
    X=0.75
    Y=0.2499

    def __init__(self, X=0, Y=0):
        self.R_HII = 0.0
        self.R_HeII = 0.0
        self.R_HeIII = 0.0
        self.n_H = 0.0
        self.n_He = 0.0
        self.M_He = 6.6464764e-24
        self.M_H = 1.6735575e-24
        self.X=X
        self.Y=Y

    #@staticmethod
    def RHS_calc(self, T, khi, g_factor):
        m_e = consts.m_e.cgs.value
        k_B = consts.k_B.cgs.value
        h = consts.h.cgs.value
        RHS = (2 * pi * m_e * k_B * T) ** 1.5 / h ** 3 * g_factor * exp(-khi / k_B / T)
        return RHS

    
    def Calc_R_HII(self, Temp):
        khi = 13.54 * u.eV
        self.R_HII = self.RHS_calc(Temp, khi.to_value(u.erg), 1 )
        return None
    
    
    def Calc_R_HeII(self, Temp):
        khi = 24.48 * u.eV
        self.R_HeII = self.RHS_calc(Temp, khi.to_value(u.erg), 4 )
        return None

    
    def Calc_R_HeIII(self, Temp):
        khi = 54.17 * u.eV
        self.R_HeIII = self.RHS_calc(Temp, khi.to_value(u.erg), 1 )
        return None

    
    def Calc_n_H(self, rho):
        self.n_H = self.X * rho / self.M_H
        return None
    
    
    def Calc_n_He(self, rho):
        self.n_He = self.Y * rho / self.M_He
        return None

    
    def Calc_consts(self,rho,Temp):
        self.Calc_n_H(rho)
        self.Calc_n_He(rho)
        self.Calc_R_HII(Temp)
        self.Calc_R_HeII(Temp)
        self.Calc_R_HeIII(Temp)
        return None


class IonizationData(tcdata.BaseData):
    """
    Class which communicates with tcdata objects.
    """
    def __init__(self):
        ionization_columnnames = ('HII_fraction','HeII_fraction','HeIII_fraction')
        super().__init__(dict(zip(ionization_columnnames,[0,0,0])),ionization_columnnames)

    @classmethod
    def initWithData(cls,x,y,z):
        obj = IonizationData()
        obj.datablock['HII_fraction'] = x
        obj.datablock['HeII_fraction'] = y
        obj.datablock['HeIII_fraction'] = z
        return obj

    def injectToDataPoint(self, datapoint_obj : tcdata.DataPoint):
        #print(self.datablock)
        
        datapoint_obj.insertColumn(self.datablock,self.column_names)
        #print(datapoint_obj.HII_fraction)
        return datapoint_obj

class Ionization:
    """
    Class mainly responsible for ionization calculations.
    """
    def __init__(self,tc_cell_object : tcdata.DataPoint,X,Y,rho=None,T=None ):
        self.Constants = ConstantsCalculator(X,Y)
        self.x = 0
        self.y = 0
        self.z = 0
        self.tc_object = tc_cell_object
        self.Constants.X = X
        self.Constants.Y = Y
        if rho is None or T is None:
            self.Constants.Calc_consts(1/self.tc_object.spec_vol,self.tc_object.temperature)      
        else:
            self.Constants.Calc_consts(rho,T)

    @staticmethod
    def diff(x1,y1,z1,x2,y2,z2):
        dx=x2-x1
        dy=y2-y1
        dz=z2-z1
        difference=sqrt(dx**2+dy**2+dz**2)
        return difference

    @staticmethod
    def __Discriminant(b,c):
        return b ** 2 - 4 *c

    def __CalcSecondOrder(self,b,c):
        det = self.__Discriminant(b,c)
        return (- b + sqrt(det))/2

    def __FirstIteration(self):
        """
        In this step we calcute x whitout y and z, calculate y without z and calculate z with y~1-z
        """
        x = y = z = 0

        b_x = self.Constants.R_HII / self.Constants.n_H
        c_x = - self.Constants.R_HII / self.Constants.n_H
        x = self.__CalcSecondOrder(b_x, c_x)

        b_y = (x* self.Constants.n_H + self.Constants.R_HeII) / self.Constants.n_He
        c_y = - self.Constants.R_HeII / self.Constants.n_He
        y = self.__CalcSecondOrder(b_y,c_y)

        b_z = (x* self.Constants.n_H + self.Constants.n_He + self.Constants.R_HeIII) / self.Constants.n_He
        c_z = - self.Constants.R_HeIII / self.Constants.n_He

        z = self.__CalcSecondOrder(b_z, c_z)
        #print(x,y,z)
        return x,y,z

    def __SecondIteration(self,x0,y0,z0):
        """
        In this step we uses previous x0 y0 and z0 values, with full n_e.
        Here we calcute frist z and after that y, because 1-y-z < 0 or should use only 1-y?
        """

        b_x = (y0 * self.Constants.n_He + 2 * z0 * self.Constants.n_He + self.Constants.R_HII) / self.Constants.n_H
        c_x = - self.Constants.R_HII / self.Constants.n_H



        b_y = (x0* self.Constants.n_H + self.Constants.R_HeII + 2 * z0 * self.Constants.n_He) / self.Constants.n_He
        c_y = (z0 - 1) * self.Constants.R_HeII / self.Constants.n_He

        y = self.__CalcSecondOrder(b_y,c_y)

        b_z = (x0* self.Constants.n_H + y * self.Constants.n_He) / (2 * self.Constants.n_He)
        c_z = - y * self.Constants.R_HeIII / (2* self.Constants.n_He)

        x = self.__CalcSecondOrder(b_x, c_x)
        
        z = self.__CalcSecondOrder(b_z, c_z)
        self.x = x
        self.y = y
        self.z = z
        return x,y,z

    def Calculation(self):
        """
        The calculation of the ionization fractions.
        """
        x,y,z = self.__FirstIteration()
        x,y,z = self.__SecondIteration(x,y,z)

        next_iterations = NewtonianIterator(x,y,z,self.Constants)
        x,y,z = next_iterations.NewtonianIteration()

        self.x = x
        self.y = y
        self.z = z

        return x,y,z

    def Reload(self, tc_cell_object : tcdata.DataPoint,X,Y,rho=None, T=None):
        self.x = 0
        self.y = 0
        self.z = 0
        self.tc_object = tc_cell_object
        self.Constants.X = X
        self.Constants.Y = Y
        if rho is None or T is None:
            self.Constants.Calc_consts(1/self.tc_object.spec_vol,self.tc_object.temperature)      
        else:
            self.Constants.Calc_consts(rho,T)

    @staticmethod
    def IonizationForRawprofile(raw_obj : tcdata.RawProfiles,X,Y):
        calculator = Ionization(raw_obj[1],X,Y)
        for cell_obj in raw_obj:
            #print(cell_obj.zone)
            if cell_obj.temperature <= 0:
                iondata = IonizationData.initWithData(0,0,0)
                iondata.injectToDataPoint(cell_obj)
                continue        
            calculator.Reload(cell_obj,X,Y)
            x,y,z = calculator.Calculation()
            iondata = IonizationData.initWithData(x,y,z)
            iondata.injectToDataPoint(cell_obj)
           # print(iondata.datablock)
           # print(cell_obj.data('HeII_fraction'))
        #print(raw_obj[15].HII_fraction)
        return raw_obj
class NewtonianIterator:
    """
    Helper class for Newton iterations.
    """

    def __init__(self,x0,y0,z0,_Constants : ConstantsCalculator):
        self.Constants = _Constants
        self.__x_vec_n = np.array([x0,y0,z0])
        self.__x_vec_np1 = np.array([x0,x0,x0])


    def __f_vector(self,x_vec) -> np.array :
        x = x_vec[0]
        y = x_vec[1]
        z = x_vec[2]

        b_x = (y * self.Constants.n_He + 2 * z * self.Constants.n_He + self.Constants.R_HII) / self.Constants.n_H
        c_x = - self.Constants.R_HII / self.Constants.n_H

        b_y = (x* self.Constants.n_H + self.Constants.R_HeII + 2 * z * self.Constants.n_He) / self.Constants.n_He
        c_y = (z - 1) * self.Constants.R_HeII / self.Constants.n_He
        b_z = (x* self.Constants.n_H + y * self.Constants.n_He) / (2 * self.Constants.n_He)
        c_z = - y * self.Constants.R_HeIII / (2* self.Constants.n_He)

        fvec = np.array(
            [x ** 2 + x * b_x  + c_x,
            y ** 2 + y * b_y + c_y,
            z ** 2 + z * b_z + c_z]
        )
        return fvec

    def __jacobian(self,x_vec = None) -> np.array :
        if x_vec is None:
            x_vec = self.__x_vec_n
        
        x = x_vec[0]
        y = x_vec[1]
        z = x_vec[2]

        b_x = (y * self.Constants.n_He + 2 * z * self.Constants.n_He + self.Constants.R_HII) / self.Constants.n_H

        b_y = (x* self.Constants.n_H + self.Constants.R_HeII + 2 * z * self.Constants.n_He) / self.Constants.n_He
        
        b_z = (x* self.Constants.n_H + y * self.Constants.n_He) / (2 * self.Constants.n_He)        

        thejacobian = np.array(
            [
            [2 *x + b_x,  x * self.Constants.n_He/self.Constants.n_H, 2* x *self.Constants.n_He / self.Constants.n_H],
            [y * self.Constants.n_H / self.Constants.n_He, 2*y + b_y, 2 * y + self.Constants.R_HeII / self.Constants.n_He] ,
            [0.5 *z * self.Constants.n_H / self.Constants.n_He, 0.5 * (z - self.Constants.R_HeIII / self.Constants.n_He), 2*z + b_z]
            ]
        )
        return thejacobian


    def NewtonianIteration(self,x_vec = None):
        """
        Starting from the third step, we use Newton-Raphson method using previous results as starting values.
        """
        if x_vec is None:
            x_vec = self.__x_vec_n
            x_vec_np1 = self.__x_vec_np1
        else:
            x_vec=np.array([x0,y0,z0])
        iteration_cnt = 0
        while True:
            Jac=self.__jacobian(x_vec)
            x_vec_np1 = x_vec - np.linalg.inv(Jac).dot(self.__f_vector(x_vec))
            
            iteration_cnt += 1
            if Ionization.diff(*x_vec_np1,*x_vec) < 1e-8 or iteration_cnt >= 1000:
                x_vec=x_vec_np1
                if iteration_cnt >= 1000:
                    raise IterationError("Iteration not converge")
                break
            x_vec=x_vec_np1
        return x_vec[0],x_vec[1],x_vec[2]

if __name__ == '__main__':


    infile = open("adat_ready.txt", "r")

    lines = [line.strip("\n") for line in infile]

    infile.close()

    data=list()

    for line in lines:
        data.append(list(float(words) for words in line.split()))


    datablock = [[] for i in range(4) ]

    print(len(lines))
    line= None

    outfile = open("kimenet.txt","w")

    # adatsor
    for line in data:
        if len(line) == 0:
            continue
        rho = 1/line[8] #adat_ready.txt
        
        T = line[10] # adat_ready.txt
        if T <= 0:
            continue
        

        datablock[0].append(T)
        # x0,y0,z0 = FirstIteration()
        # x0,y0,z0 = SecondIteration(x0,y0,z0)
        # x0,y0,z0 = NewtonianIteration(x0,y0,z0)
        ionization_obj = Ionization(tcdata.DataPoint(""),0.75,0.2499,rho,T)
        x0,y0,z0 = ionization_obj.Calculation()
        datablock[1].append(x0)
        datablock[2].append(y0)
        datablock[3].append(z0)
        
        outfile.write("{0} {1:8.6E} {2:8.6E} {3:8.6E} {4:8.6E} {5:8.6E} {6:8.6E}\n".format(T,x0,y0,z0,line[2],line[3],line[4]))

        ###Checking

        

        """n_e=x2*Constants.n_H + (y2 + 2 * z2) * Constants.n_He

        n_e0=x0*Constants.n_H + (y0 + 2 * z0) * Constants.n_He

        egy1 = x2 * n_e / (1 - x2) - Constants.R_HII
        egy2 = y2 * n_e / (1 - y2 - z2) - Constants.R_HeII
        egy3 = z2 * n_e / y2 - Constants.R_HeIII

        egy01 = x0 * n_e0 / (1 - x0) - Constants.R_HII
        egy02 = y0 * n_e0 / (1 - y0 - z0) - Constants.R_HeII
        egy03 = z0 * n_e0 / y0 - Constants.R_HeIII
        """
        #print("            pontosság: ", (egy1/Constants.R_HII+egy1/(x0 * n_e / (1 - x0)))/2)

        #y0 = #fsolve(f02, 0.675486)
        #z0 = #fsolve(f03, 0.664258)
        #print(x0,y0,z0)
        #x,y,z, = -1, -1,-1
        #if T > 100000:
        #    x0,y0,z0 = 0.99 , 0.01, 0.99
        #elif T > 50000:
        #    x0,y0,z0 = 0.9576, 0.6548, 0.4012
        #elif T > 20000:
        #    x0,y0,z0 = 0.9, 0.25, 0.01
        #elif T > 10000:
        #    x0, y0, z0 = 0.8, 0.01, 0.001
        #else:
        #    x0,y0,z0 = 0.2, 0.01, 0.001

        res = []
        #while x < 0 or x > 1 or y < 0 or y > 1 or z < 0 or z > 1:
            #x,y,z = fsolve(func, [x0,y0,z0])
            #print(x,y,z)
        #for x0 in [ j * 0.1 for j in range(10)]:
        #    for y0 in [ j * 0.1 for j in range(10)]:
        #        for z0 in [ j * 0.1 for j in range(10)]:
        #            x,y,z = fsolve(func, [x0,y0,z0])
        #            if x >= 0 and x <= 1 and y >= 0 and y <= 1 and z >= 0 and z <= 1:
        #                res.append([x,y,z])
        #                print(x,y,z)


        #print(res)    
        #break
        #print(rho,T,x0)#,y0,z0)
    outfile.close()

    print(type(datablock[1][0]))

    #exit()
    from matplotlib import pyplot as plt

    # plt.plot(datablock[0],datablock[1], label = "x")
    # plt.plot(datablock[0],datablock[2], label = 'y')
    # plt.plot(datablock[0],datablock[3], label = 'z')
    # plt.legend()
    # plt.xlim(6000,1e5)
    # plt.show()

    a_point = tcdata.DataPoint("45\t136\t12\t50\t5.6\t2e5\t4.13e+29\t125.68\t1.25e40\t1.2e38\t12000\t1e29\t12e5\t1.5e4\t1e2\t0.98\n")
    ions_from_tc = Ionization(a_point,0.75,0.2496,1e-2,1e5)
    x,y,z =ions_from_tc.Calculation()
    iondata = IonizationData.initWithData(x,y,z)
    a_point = iondata.injectToDataPoint(a_point)
    print(a_point.datablock)
    fort19_path='/home/gabesz/SPHERLS_playground/bpf_konv/konv-b/fort.19'
    fort19_data=tcdata.RawProfiles(fort19_path)
    fort19_data = Ionization.IonizationForRawprofile(fort19_data,0.75,0.2496)
    print(fort19_data[15].spec_vol)
    fort19_handler = tcdata.LimitCycle(fort19_data)

    for i in range(len(fort19_handler.profiles[2].HeII_fraction)):
        print(fort19_handler.profiles[2].zone[i],fort19_handler.profiles[2].temperature[i],fort19_handler.profiles[2].HeII_fraction[i])

    plt.plot(fort19_handler.profiles[2].temperature[2:147],fort19_handler.profiles[2].HII_fraction[2:147], label= 'HII')
    plt.plot(fort19_handler.profiles[2].temperature[2:147],fort19_handler.profiles[2].HeII_fraction[2:147], label= 'HeII')
    plt.plot(fort19_handler.profiles[2].temperature[2:147],fort19_handler.profiles[2].HeIII_fraction[2:147], label= 'HeIII')
    plt.legend()
    plt.xlim(6500,1e5)
    plt.show()
