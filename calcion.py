#Calculate ionisation fractions
#
# v 0.1
#TODO:
# - Refactor iterations to OOP
# - interface with tcdata module
# - Callable functions from outside
# - Testing utilites
# - Deleting unneded imports and ineffective attempts
# - Refactor again with the working tests above


from math import pi
from math import exp
from math import sqrt
import astropy.constants as consts
import astropy.units as u
from astropy.units.core import _recreate_irreducible_unit
import numpy as np

class IterationError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class Constants:
    R_HII = 0.0
    R_HeII = 0.0
    R_HeIII = 0.0
    n_H = 0.0
    n_He = 0.0
    M_He = 6.6464764e-24
    M_H = 1.6735575e-24
    X=0.75
    Y=0.2499

    #@staticmethod
    def RHS_calc(self, T, khi, g_factor):
        m_e = consts.m_e.cgs.value
        k_B = consts.k_B.cgs.value
        h = consts.h.cgs.value
        RHS = (2 * pi * m_e * k_B * T) ** 1.5 / h ** 3 * g_factor * exp(-khi / k_B / T)
        return RHS

    @classmethod
    def Calc_R_HII(cls, Temp):
        khi = 13.54 * u.eV
        Constants.R_HII = cls.RHS_calc(cls , Temp, khi.to_value(u.erg), 1 )
        return None
    
    @classmethod
    def Calc_R_HeII(cls, Temp):
        khi = 24.48 * u.eV
        Constants.R_HeII = cls.RHS_calc(cls, Temp, khi.to_value(u.erg), 4 )
        return None

    @classmethod
    def Calc_R_HeIII(cls, Temp):
        khi = 54.17 * u.eV
        Constants.R_HeIII = cls.RHS_calc(cls, Temp, khi.to_value(u.erg), 1 )
        return None

    @classmethod
    def Calc_n_H(cls, rho):
        cls.n_H = cls.X * rho / cls.M_H
        return None
    
    @classmethod
    def Calc_n_He(cls, rho):
        cls.n_He = cls.Y * rho / cls.M_He
        return None

    @classmethod
    def Calc_consts(cls,rho,Temp):
        cls.Calc_n_H(rho)
        cls.Calc_n_He(rho)
        cls.Calc_R_HII(Temp)
        cls.Calc_R_HeII(Temp)
        cls.Calc_R_HeIII(Temp)
        return None


def diff(x1,y1,z1,x2,y2,z2):
    dx=x2-x1
    dy=y2-y1
    dz=z2-z1
    difference=sqrt(dx**2+dy**2+dz**2)
    return difference



Constants.Calc_R_HeII(15000)

infile = open("adat_ready.txt", "r")

lines = [line.strip("\n") for line in infile]

infile.close()

data=list()

for line in lines:
    data.append(list(float(words) for words in line.split()))

def Discriminant(b,c):
    return b ** 2 - 4 *c

def CalcSecondOrder(b,c):
    det = Discriminant(b,c)
    return (- b + sqrt(det))/2

def FirstIteration():
    """
    In this step we calcute x whitout y and z, calculate y without z and calculate z with y~1-z
    """
    x = y = z = 0

    b_x = Constants.R_HII / Constants.n_H
    c_x = - Constants.R_HII / Constants.n_H
    x = CalcSecondOrder(b_x, c_x)

    b_y = (x* Constants.n_H + Constants.R_HeII) / Constants.n_He
    c_y = - Constants.R_HeII / Constants.n_He
    y = CalcSecondOrder(b_y,c_y)

    b_z = (x* Constants.n_H + Constants.n_He + Constants.R_HeIII) / Constants.n_He
    c_z = - Constants.R_HeIII / Constants.n_He

    z = CalcSecondOrder(b_z, c_z)
    #print(x,y,z)
    return x,y,z

def SecondIteration(x0,y0,z0):
    """
    In this step we uses previous x0 y0 and z0 values, with full n_e.
    Here we calcute frist z and after that y, because 1-y-z < 0 or should use only 1-y?
    """

    b_x = (y0 * Constants.n_He + 2 * z0 * Constants.n_He + Constants.R_HII) / Constants.n_H
    c_x = - Constants.R_HII / Constants.n_H



    b_y = (x0* Constants.n_H + Constants.R_HeII + 2 * z0 * Constants.n_He) / Constants.n_He
    c_y = (z0 - 1) * Constants.R_HeII / Constants.n_He

    y = CalcSecondOrder(b_y,c_y)

    b_z = (x0* Constants.n_H + y * Constants.n_He) / (2 * Constants.n_He)
    c_z = - y * Constants.R_HeIII / (2* Constants.n_He)

    x = CalcSecondOrder(b_x, c_x)
    
    z = CalcSecondOrder(b_z, c_z)

    return x,y,z

def f_vector(x_vec) -> np.array :
    x = x_vec[0]
    y = x_vec[1]
    z = x_vec[2]

    b_x = (y * Constants.n_He + 2 * z * Constants.n_He + Constants.R_HII) / Constants.n_H
    c_x = - Constants.R_HII / Constants.n_H

    b_y = (x* Constants.n_H + Constants.R_HeII + 2 * z * Constants.n_He) / Constants.n_He
    c_y = (z - 1) * Constants.R_HeII / Constants.n_He
    b_z = (x* Constants.n_H + y * Constants.n_He) / (2 * Constants.n_He)
    c_z = - y * Constants.R_HeIII / (2* Constants.n_He)

    fvec = np.array(
        [x ** 2 + x * b_x  + c_x,
        y ** 2 + y * b_y + c_y,
        z ** 2 + z * b_z + c_z]
    )
    return fvec

def jacobian(x_vec) -> np.array :
    x = x_vec[0]
    y = x_vec[1]
    z = x_vec[2]

    b_x = (y * Constants.n_He + 2 * z * Constants.n_He + Constants.R_HII) / Constants.n_H

    b_y = (x* Constants.n_H + Constants.R_HeII + 2 * z * Constants.n_He) / Constants.n_He
    
    b_z = (x* Constants.n_H + y * Constants.n_He) / (2 * Constants.n_He)
    

    thejacobian = np.array(
        [
          [2 *x + b_x,  x * Constants.n_He/Constants.n_H, 2* x *Constants.n_He / Constants.n_H],
          [y * Constants.n_H / Constants.n_He, 2*y + b_y, 2 * y + Constants.R_HeII / Constants.n_He] ,
          [0.5 *z * Constants.n_H / Constants.n_He, 0.5 * (z - Constants.R_HeIII / Constants.n_He), 2*z + b_z]
        ]
    )
    return thejacobian


def NextIteration(x0,y0,z0):
    """
    Here are the iterations supposed to use more than once.
    """
    b_x = (y0 * Constants.n_He + 2 * z0 * Constants.n_He + Constants.R_HII) / Constants.n_H
    c_x = - Constants.R_HII / Constants.n_H

    b_y = (x0* Constants.n_H + Constants.R_HeII + 2 * z0 * Constants.n_He) / Constants.n_He
    c_y = (z0 - 1) * Constants.R_HeII / Constants.n_He
    b_z = (x0* Constants.n_H + y0 * Constants.n_He) / (2 * Constants.n_He)
    c_z = - y0 * Constants.R_HeIII / (2* Constants.n_He)

    x = CalcSecondOrder(b_x, c_x)
    y = CalcSecondOrder(b_y,c_y)
    z = CalcSecondOrder(b_z, c_z)
    return x,y,z

#print(data[15])

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
    Constants.Calc_consts(rho, T)

    datablock[0].append(T)
    x0,y0,z0 = FirstIteration()
    x0,y0,z0 = SecondIteration(x0,y0,z0)
    x_vec=np.array([x0,y0,z0])
    iteration_cnt = 0
    while True:
        Jac=jacobian(x_vec)
        x_vec_new = x_vec - np.linalg.inv(Jac).dot(f_vector(x_vec))
        print(diff(*x_vec_new, *x_vec))
        iteration_cnt += 1
        if diff(*x_vec_new,*x_vec) < 1e-8 or iteration_cnt >= 1000:
            x_vec=x_vec_new
            print(iteration_cnt)
            if iteration_cnt >= 1000:
                raise IterationError("Iteration not converge")
            break
        x_vec=x_vec_new
    x0=x_vec[0]
    y0=x_vec[1]
    z0=x_vec[2]
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

plt.plot(datablock[0],datablock[1], label = "x")
plt.plot(datablock[0],datablock[2], label = 'y')
plt.plot(datablock[0],datablock[3], label = 'z')
plt.legend()
plt.xlim(6000,1e5)
plt.show()