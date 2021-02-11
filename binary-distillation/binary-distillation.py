import numpy as np
import matplotlib.pyplot as plt

#*******************************************************************************
#This program calculates the boiling point versus gas and liquid phase
#compositions for nonideal binary mixture. It takes as input the Antoine constants of
#both the components and the temperature and composition of the azeotrope.
#The activities coefficients are calculated through the Van Laar equations
#*******************************************************************************

#Antoine equation: https://en.wikipedia.org/wiki/Antoine_equation
def antoineEq(A, B, C, t):
    P_0i=A -B/(t+C)
    P_0i=10**P_0i
    return P_0i

#Determination of the boiling point using the Antoine quation
def boilingPoint(P,A,B,C):
     t= B/(A-np.log10(P))-C
     return t

#Determination of the activity coefficients at the azeotropic point
def azeotropeActivityCoefficient(P, P0):
    gamma= P/P0
    return gamma

#Determination of the constants A and B in Van Laar equations
def VanLaarConstantsAB(X1_az,X2_az,gamma1,gamma2):
    A_laar= np.log(gamma1)*(1+(X2_az*np.log(gamma2)/(X1_az*np.log(gamma1))))**2
    B_laar= np.log(gamma2)*(1+(X1_az*np.log(gamma1)/(X2_az*np.log(gamma2))))**2
    return A_laar,B_laar

#Determination of the activity coefficients using the Van Laar Equations
def VanLaarActivityCoefficients(A_laar,B_laar,X1,X2):
    gammaL1=A_laar/(1+X1/X2*A_laar/B_laar)**2
    gammaL1= np.exp(gammaL1)
    gammaL2=B_laar/(1+X2/X1*B_laar/A_laar)**2
    gammaL2= np.exp(gammaL2)
    return gammaL1,gammaL2

#Dalton's law of partial pressures.
def DaltonLaw(Ptot,gammaL1, gammaL2,X1,X2,P01,P02):
    Psys=gammaL1*X1*P01+gammaL2*X2*P02
    return Ptot-Psys

#*******************************************************************************
#A1,B1,C1 = Antoine Constants for component 1
#A2,B2,C2 = Antoine Constants for component 2
#azeotropeT= azeotrope temperature (Celsius)
#X1_az,X2_az = azeotrope composition (Molar fractions: X1,X2)
#Ptot= atmospheric pressure (torr)
#P_01azeotrope, P_02azeotrope= partial vapor pressures of component 1 and
#component 2 at the azeotropic point.
#gamma1_azeotrope, gamma2_azeotrope= activity coefficients of component 1
# and component 2 at the azeotropic point
#A_vanLaar,B_vanLaar= Van Laar coefficients
#BP_comp1, BP_comp2 = boiling points of component 1 and component 2
#*******************************************************************************
ComponentA=input("Enter the name of the component A: ")
A1,B1,C1 = [float(x) for x in input("Enter the Antoine Constants for component 1 (A,B,C): ").split()]
ComponentB=input("Enter the name of the component B: ")
A2,B2,C2 = [float(x) for x in input("Enter the Antoine Constants for component 2 (A,B,C): ").split()]
azeotropeT=float(input("Enter the azeotrope temperature (Celsius): "))
X1_az,X2_az = [float(x) for x in input("Enter azeotrope composition (X1,X2): ").split()]

# The calculations refer to the pressure of 1 atm (760 torr)
Ptot=760.0

P_01azeotrope= antoineEq(A1, B1, C1, azeotropeT)
P_02azeotrope= antoineEq(A2, B2, C2, azeotropeT)
gamma1_azeotrope=azeotropeActivityCoefficient(Ptot,P_01azeotrope)
gamma2_azeotrope=azeotropeActivityCoefficient(Ptot,P_02azeotrope)
A_vanLaar,B_vanLaar=VanLaarConstantsAB(X1_az,X2_az,gamma1_azeotrope,gamma2_azeotrope)
BP_comp1 =boilingPoint(Ptot,A1,B1,C1)
BP_comp2 =boilingPoint(Ptot,A2,B2,C2)

#set y axis limits
if BP_comp1>BP_comp2:
    lim_yup=BP_comp1+2.5
else:
    lim_yup=BP_comp2+2.5
print (BP_comp1,BP_comp2)
if azeotropeT<BP_comp1 and azeotropeT<BP_comp2:
    lim_ydown=azeotropeT-1
elif BP_comp1>BP_comp2:
    lim_yup=azeotropeT+2.5
    lim_ydown=BP_comp2-1
elif BP_comp1<BP_comp2:
    lim_yup=azeotropeT+2.5
    lim_ydown=BP_comp1-1


#Estimation of the upper temperature limit for the bisection method
if azeotropeT<BP_comp1 and azeotropeT<BP_comp2:
   if BP_comp1>BP_comp2:
       tlim=BP_comp1+BP_comp1/2
   else:
       tlim=BP_comp2+BP_comp2/2
else:
    tlim=azeotropeT+azeotropeT/2

#Error tolerance in Dalton law
epsilon=0.001
#Max number of bisection iterations
maxiter= 1000
#Initialization of the temperature array
Tcurve=np.zeros(101)
# Molar fractions componet 1
X= np.arange(0.00,1.01,0.01)
# Molar fractions componet 2
X_2= np.arange(1.00,-0.01,-0.01)
# Initialization of the array for the molar fractions of a component in vapor phase
XP=np.zeros(101)

for i in range (0, len(X)):

    if X[i]==0 and X_2[i]==1:
        Tcurve[i]=BP_comp2
        XP[i]=0

    elif X[i]==1 and X_2[i]==0:
        Tcurve[i]=BP_comp1
        XP[i]=1
    else:
        gammaL1,gammaL2= VanLaarActivityCoefficients(A_vanLaar,B_vanLaar,X[i],X_2[i])

# Bisection method to find the temperature that satisfies Dalton Law
        t0=0
        te=tlim
        for j in range (0,maxiter):
            t= (t0+te)/2
            P1= antoineEq(A1, B1, C1, t)
            P2= antoineEq(A2, B2, C2, t)
            Delta=DaltonLaw(Ptot,gammaL1, gammaL2,X[i],X_2[i],P1,P2)
            if abs(Delta)< epsilon:
               Tcurve[i]=t
               XP[i]=P1*gammaL1*X[i]/(P1*gammaL1*X[i]+P2*gammaL2*X_2[i])
               break
            elif Delta<0:
               te=t
            else:
               t0=t
#*******************************************************************************
#Plot T vs XA(molar fraction of component 1)
#*******************************************************************************
plt.title("Vapor-liquid equilibrium (1atm): " + ComponentA + " - "+ ComponentB)
plt.plot( X,Tcurve, '#00004f',label='XA liquid')
plt.plot( XP,Tcurve, 'red',label='XA vapor')
plt.xlabel ("XA - "+ ComponentA)
plt.ylabel ("Temperature (Degrees Celsius)")
plt.ylim([lim_ydown,lim_yup])
plt.xlim([0,1])
plt.fill_between(XP, Tcurve, lim_yup, color='orange', label='vapor phase')
plt.fill_between(XP, Tcurve,color='pink',label='vapor-liquid equilibrium')
plt.fill_between(X, Tcurve,color='lightblue',label='liquid phase')
plt.legend(loc="upper left")
if azeotropeT<BP_comp1 and azeotropeT<BP_comp2:
    plt.annotate('azeotropic point', xy=(X1_az,azeotropeT),xytext=(X1_az, azeotropeT+2),horizontalalignment="center",arrowprops={'arrowstyle': '->'})
else:
    plt.annotate('azeotropic point', xy=(X1_az,azeotropeT),xytext=(X1_az, azeotropeT-2),horizontalalignment="center",arrowprops={'arrowstyle': '->'})
plt.show()

#*******************************************************************************
#Plot T vs XB(molar fraction of component 2)
#*******************************************************************************
for i in range(0,len(XP)):
   XP[i]=1-XP[i]
plt.title("Vapor-liquid equilibrium (1atm): "+ ComponentA + " - "+ ComponentB)
plt.plot( X_2,Tcurve, '#00004f',label='XB liquid')
plt.plot( XP,Tcurve, 'red',label='XB vapor')
plt.xlabel ("XB - "+ ComponentB)
plt.ylabel ("Temperature (Degrees Celsius)")
plt.ylim([lim_ydown,lim_yup])
plt.xlim([0,1])
plt.fill_between(XP, Tcurve, lim_yup, color='orange', label='vapor phase')
plt.fill_between(XP, Tcurve,color='pink',label='vapor-liquid equilibrium')
plt.fill_between(X_2, Tcurve,color='lightblue', label='liquid phase')
plt.legend(loc="upper left")
if azeotropeT<BP_comp1 and azeotropeT<BP_comp2:
    plt.annotate('azeotropic point', xy=(X2_az,azeotropeT),xytext=(X2_az, azeotropeT+2),horizontalalignment="center",arrowprops={'arrowstyle': '->'})
else:
    plt.annotate('azeotropic point', xy=(X2_az,azeotropeT),xytext=(X2_az, azeotropeT-2),horizontalalignment="center",arrowprops={'arrowstyle': '->'})
plt.show()

#*******************************************************************************
#Examples of input:
#
#Enter the name of the component A: Ethyl Acetate
#Enter the Antoine Constants for component 1 (A,B,C): 7.10179 1244.951 217.881
#Enter the name of the component B: Cyclohexane
#Enter the Antoine Constants for component 2 (A,B,C): 6.84941 1206.001 223.148
#Enter the azeotrope temperature (Celsius): 72.8
#Enter azeotrope composition (X1,X2): 0.531 0.469

#Enter the name of the component A: Hexane
#Enter the Antoine Constants for component 1 (A,B,C): 6.88555 1175.817 224.867
#Enter the name of the component B:  Ethyl Acetate
#Enter the Antoine Constants for component 2 (A,B,C): 7.10179 1244.951 217.881
#Enter the azeotrope temperature (Celsius): 65.2
#Enter azeotrope composition (X1,X2): 0.606 0.394

#Enter the name of the component A: 1-Propano
#Enter the Antoine Constants for component 1 (A,B,C): 7.74416 1437.686 198.463
#Enter the name of the component B: Cyclohexane
#Enter the Antoine Constants for component 2 (A,B,C): 6.84941 1206.001 223.148
#Enter the azeotrope temperature (Celsius): 74.7
#Enter azeotrope composition (X1,X2): 0.241 0.759
