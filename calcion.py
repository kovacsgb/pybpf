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

from scipy.optimize import fsolve
import sympy as sympy
from math import pi
from math import exp
from math import sqrt
import astropy.constants as consts
import astropy.units as u

x, y, z = sympy.symbols("x y z")
c1,c2,c3,c4,c5 = sympy.symbols("c1 c2 c3 c4 c5")

#print(consts.h.cgs.value,consts.k_B.cgs.value)

#print((2*pi * consts.m_e.cgs.value * consts.k_B.cgs.value * 1000) ** 1.5 / (consts.h.cgs.value ** 3))


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


def Iteration(x1,y1,z1):
    #print("It start:",x1,y1,z1)
    z2 = 0
    y2 = 0
    x2 = 0

    for i in range(100):
        bval=((y1 + 2 *z1) * Constants.n_He +Constants.R_HII)/Constants.n_H
        cval= - Constants.R_HII / Constants.n_H
        det = bval ** 2 - 4 * cval
        x2 = ( - (1 * bval) + sqrt(det)) / 2


        bval_z = (Constants.n_H * x1 + Constants.n_He + Constants.R_HeIII) / ( Constants.n_He)
        cval_z = - ( Constants.R_HeIII)  / ( Constants.n_He)
        
        det3 = (bval_z ** 2 - 4 * cval_z )
        z2=  ( - (1 * bval_z) + sqrt(det3)) / 2

        bval_y=(Constants.n_H * x1 + Constants.R_HeII + 2 * z1 * Constants.n_He) / Constants.n_He
        cval_y =(z1 * Constants.R_HeII - Constants.R_HeII ) / Constants.n_He
        det2 = (bval_y ** 2 - 4 * cval_y )
        y2=  ( - (1 * bval_y) + sqrt(det2)) / 2

        #print("It:", i, x2,y2,z2)
        

        ###step 2:
        bval=((y2 + 2 *z2) * Constants.n_He +Constants.R_HII)/Constants.n_H
        cval= - Constants.R_HII / Constants.n_H
        det = bval ** 2 - 4 * cval
        x1 = ( - (1 * bval) + sqrt(det)) / 2

        #bval_z = (Constants.n_H * x2 + y2 * Constants.n_He ) / (2* Constants.n_He)
        #cval_z = - (y2 * Constants.R_HeIII)  / (2* Constants.n_He)
        bval_z = (Constants.n_H * x2 + Constants.n_He + Constants.R_HeIII) / ( Constants.n_He)
        cval_z = - ( Constants.R_HeIII)  / ( Constants.n_He     )   
        
        det3 = (bval_z ** 2 - 4 * cval_z )
        z1=  ( - (1 * bval_z) + sqrt(det3)) / 2

        bval_y=(Constants.n_H * x2 + Constants.R_HeII + 2 * z2 * Constants.n_He) / Constants.n_He
        cval_y =(z2 * Constants.R_HeII - Constants.R_HeII ) / Constants.n_He
        det2 = (bval_y ** 2 - 4 * cval_y )
        y1=  ( - (1 * bval_y) + sqrt(det2)) / 2


        

        #print(det3, x2, y2, z2)


    return x1, y1, z1




Constants.Calc_R_HeII(15000)


def f1(x,y,z):
    result = x * (Constants.n_H * x + Constants.n_He * y + Constants.n_He * 2 * z) / (1 - x) - Constants.R_HII
    return result

def f2(x,y,z):
    result = y * (Constants.n_H * x + Constants.n_He * y + Constants.n_He * 2 * z) / (1 - y - z) - Constants.R_HeII
    return result

def f3(x,y,z):
    result = z *(Constants.n_H * x + Constants.n_He * y + Constants.n_He * 2 * z) / y - Constants.R_HeIII
    return result

def f01(x):
    result = x ** 2 / (1 - x) - Constants.R_HII
    return result

def f02(y):
    result = y ** 2 / (1 - y) - Constants.R_HeII
    return result

def f03(z):
    result = z ** 2 / (1 - z) - Constants.R_HeIII
    return result



def func(x):
    res1=f1(x[0],x[1],x[2])
    res2=f2(x[0],x[1],x[2])
    res3=f3(x[0],x[1],x[2])
    return [res1, res2, res3]

#beolvasás

#infile = open("fort.19", "r")

infile = open("adat_ready.txt", "r")

lines = [line.strip("\n") for line in infile]

infile.close()

data=list()

for line in lines:
    data.append(list(float(words) for words in line.split()))




#print(data[15])

print(len(lines))
line= None

outfile = open("kimenet.txt","w")

# adatsor
for line in data:
    if len(line) == 0:
        continue
    #print(line[6]/line[17],line[5])
    #rho = line[6]/line[17]
    rho = 1/line[8] #adat_ready.txt
    
    #T = line[5]
    T = line[10] # adat_ready.txt
    #print(line)
   # print(T,"{0:E} {1:E} {2:E}".format(rho, line[12], line[8]))
    
    if T <= 0:
        continue
    Constants.Calc_consts(rho, T)
    kif= x ** 2 / (1 - x) - Constants.R_HII / Constants.n_H

    bval=Constants.R_HII/Constants.n_H

    det = bval ** 2 + 4 * bval


    #Érdemes lehet felcserélni z és y számolását. (z magasabb hőmérsékleten lesz 1, és 1-z jobban közelíti ekkor y-t.) Ezzel csinálni az első iterációs lépést.
    x0 = ( - (1 * bval) + sqrt(det)) / 2
    x02 = ( - (1 * bval) - sqrt(det)) / 2
    #print(T," K and x0=",x0," ",x02, "and det =",bval ** 2 + 4 * bval, "bval =",bval)
   # print( sympy.nsolve(kif,x0))

    ###calc z
    bval_z = (Constants.n_H * x0 + Constants.n_He + Constants.R_HeIII) / ( Constants.n_He)
    cval_z = - (Constants.R_HeIII)  / ( Constants.n_He)
    det3 = (bval_z ** 2 - 4 * cval_z )
    z0=  ( - (1 * bval_z) + sqrt(det3)) / 2
    z02=  ( - (1 * bval_z) - sqrt(det3)) / 2

    ###calc y
    bval_y = (Constants.n_H * x0 + Constants.R_HeII + 2 * z0 * Constants.n_He) / Constants.n_He
    cval_y =(z0 * Constants.R_HeII - Constants.R_HeII ) / Constants.n_He
    det2 = (bval_y ** 2 - 4 * cval_y )
    y0=  ( - (1 * bval_y) + sqrt(det2)) / 2
    y02=  ( - (1 * bval_y) - sqrt(det2)) / 2
    #print("          and y0=",y0," ",y02, "and det2 =",det2, "bval_y=",bval_y, "cval_y=", cval_y)

   # print("          and z0=",z0," ",z02, "and det3 =",det3, "bval_z=",bval_z, "cval_z=", cval_z)


    ###Second iteration
    bval=((y0 + 2 *z0) * Constants.n_He +Constants.R_HII)/Constants.n_H
    cval= - Constants.R_HII / Constants.n_H
    det = bval ** 2 - 4 * cval
    x1 = ( - (1 * bval) + sqrt(det)) / 2
    
    bval_z = (Constants.n_H * x0 + y0 * Constants.n_He ) / (2* Constants.n_He)
    cval_z = - (y0 * Constants.R_HeIII)  / (2* Constants.n_He)
    
    det3 = (bval_z ** 2 - 4 * cval_z )

    z1=  ( - (1 * bval_z) + sqrt(det3)) / 2  

    bval_y=(Constants.n_H * x0 + Constants.R_HeII + 2 * z0 * Constants.n_He) / Constants.n_He
    cval_y =(z0 * Constants.R_HeII - Constants.R_HeII ) / Constants.n_He
    det2 = (bval_y ** 2 - 4 * cval_y )
    y1=  ( - (1 * bval_y) + sqrt(det2)) / 2

  

   ### Third iteration
    bval=((y1 + 2 *z1) * Constants.n_He +Constants.R_HII)/Constants.n_H
    cval= - Constants.R_HII / Constants.n_H
    det = bval ** 2 - 4 * cval
    x2 = ( - (1 * bval) + sqrt(det)) / 2

    bval_z = (Constants.n_H * x1 + y1 * Constants.n_He ) / (2* Constants.n_He)
    cval_z = - (y1 * Constants.R_HeIII)  / (2* Constants.n_He)
    
    det3 = (bval_z ** 2 - 4 * cval_z )
    z2=  ( - (1 * bval_z) + sqrt(det3)) / 2 


    bval_y=(Constants.n_H * x1 + Constants.R_HeII + 2 * z1 * Constants.n_He) / Constants.n_He
    cval_y =(z1 * Constants.R_HeII - Constants.R_HeII ) / Constants.n_He
    det2 = (bval_y ** 2 - 4 * cval_y )
    y2=  ( - (1 * bval_y) + sqrt(det2)) / 2

 
    x3,y3,z3 = Iteration(x2,y2,z2)
    #x3 = x2
    #y3 = y2
    #z3 = z2
   ### Last Iteration
    bval=((y3 + 2 *z3) * Constants.n_He +Constants.R_HII)/Constants.n_H
    cval= - Constants.R_HII / Constants.n_H
    det = bval ** 2 - 4 * cval
    x2 = ( - (1 * bval) + sqrt(det)) / 2

    bval_z = (Constants.n_H * x3 + y3 * Constants.n_He ) / (2* Constants.n_He)
    cval_z = - (y3 * Constants.R_HeIII)  / (2* Constants.n_He)
    
    det3 = (bval_z ** 2 - 4 * cval_z )
    z2=  ( - (1 * bval_z) + sqrt(det3)) / 2 


    bval_y=(Constants.n_H * x2 + Constants.R_HeII + 2 * z2 * Constants.n_He) / Constants.n_He
    cval_y =(z2 * Constants.R_HeII - Constants.R_HeII ) / Constants.n_He
    det2 = (bval_y ** 2 - 4 * cval_y )
    y2=  ( - (1 * bval_y) + sqrt(det2)) / 2



    #print(T, " K, x = ", x2, " x0 = ", x1, " Delta X = ", (x3-x2))
    #print("            y = ", y2, " y0 = ", y1, " Delta Y = ", (y3-y2))
    #print("            z = ", z2, " z0 = ", z1, " Delta Z = ", (z3-z2) )

    #print(T, " K, x = ", x2, " y = ", y2, " z = ", z2)
    #print("{0} K, x = {1:8.6E}\ty = {2:8.6E}\tz = {3:8.6E} ;;; xe = {4:8.6E}, ye = {5:8.6E}, ze = {6:8.5E}".format(T,x2,y2,z2,line[2],line[3],line[4]))
    #print("{0} K, x = {1:8.6E}\ty = {2:8.6E}\tz = {3:8.6E} ;;; dx = {4:8.6E}, dy = {5:8.6E}, dz = {6:8.5E}".format(T,x2,y2,z2,line[2]-x2,line[3]-y2,line[4]-z2))
    #print("{0} K, x = {1:8.6E}\ty = {2:8.6E}\tz = {3:8.6E} ;;; dx = {4:8.6E}, dy = {5:8.6E}, dz = {6:8.5E}".format(T,x2,y2,z2,abs(line[2]-x2)/line[2],abs((line[3]-y2)/line[3]),abs((line[4]-z2)/line[4])))
    
    
    outfile.write("{0} {1:8.6E} {2:8.6E} {3:8.6E} {4:8.6E} {5:8.6E} {6:8.6E}\n".format(T,x2,y2,z2,line[2],line[3],line[4]))

    ###Checking

    

    n_e=x2*Constants.n_H + (y2 + 2 * z2) * Constants.n_He

    n_e0=x0*Constants.n_H + (y0 + 2 * z0) * Constants.n_He

    egy1 = x2 * n_e / (1 - x2) - Constants.R_HII
    egy2 = y2 * n_e / (1 - y2 - z2) - Constants.R_HeII
    egy3 = z2 * n_e / y2 - Constants.R_HeIII

    egy01 = x0 * n_e0 / (1 - x0) - Constants.R_HII
    egy02 = y0 * n_e0 / (1 - y0 - z0) - Constants.R_HeII
    egy03 = z0 * n_e0 / y0 - Constants.R_HeIII

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