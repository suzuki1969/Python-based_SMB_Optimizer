# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# ASMB model was constructed by Kensuke Suzuki in 2020
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

import numpy as np

# fill the number of dataset in CSV files 
NData = 30 # Number Of Dataset

###################Fru, Gruの順#######################

filename = "SubStepTime.csv"
ALLHT = np.loadtxt(filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

HT1 = ALLHT[0,0:NData]
HT2 = ALLHT[1,0:NData]
HT3 = ALLHT[2,0:NData]
HT4 = ALLHT[3,0:NData]
HT5 = ALLHT[4,0:NData]
HT6 = ALLHT[5,0:NData]
StepTimeInit = np.zeros(NData)
for i in range(NData):
    StepTimeInit[i] = sum(ALLHT[0:5,i]) 

filename = "SubStepVelocity.csv"
ALLU = np.loadtxt(filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

URexp  = ALLU[0,0:NData]
URaffexp = ALLU[1,0:NData]
UFeedRaffexp = ALLU[2,0:NData]
UExtexp = ALLU[3,0:NData]
Umin = 0.000001

filename = "ProductConcentration.csv"
ALLConc = np.loadtxt(filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

Concentration_of_Feed_A = ALLConc[0,0:NData]
Concentration_of_Feed_B = ALLConc[1,0:NData]
Concentration_of_Raffinate_A = ALLConc[2,0:NData]
Concentration_of_Raffinate_B = ALLConc[3,0:NData]
Concentration_of_Extract_A = ALLConc[4,0:NData]
Concentration_of_Extract_B = ALLConc[5,0:NData]

CFeed = {}
CExtract = {}
CRaffinate = {}

for i in range(NData):
    CFeed[i+1,1] = Concentration_of_Feed_A[i]
    CExtract[i+1,1] = Concentration_of_Extract_A[i]
    CRaffinate[i+1,1] = Concentration_of_Raffinate_A[i]
    CFeed[i+1,2] = Concentration_of_Feed_B[i]
    CExtract[i+1,2] = Concentration_of_Extract_B[i]
    CRaffinate[i+1,2] = Concentration_of_Raffinate_B[i]

# fluid velocity of Column 1
U1Init = np.array([[URexp[0],URaffexp[0],Umin,Umin,Umin,URexp[0],UExtexp[0],URexp[0],URexp[0],URexp[0],URexp[0],URexp[0]]])
# fluid velocity of Column 2
U2Init = np.array([[URexp[0],URaffexp[0],Umin,Umin,Umin,URexp[0],Umin,URexp[0],URexp[0],URexp[0],URexp[0],URexp[0]]])
# fluid velocity of Column 3
U3Init = np.array([[URexp[0],URaffexp[0],UFeedRaffexp[0],UFeedRaffexp[0],UFeedRaffexp[0],URexp[0],Umin,URexp[0],URexp[0],URexp[0],URexp[0],URexp[0]]])
# fluid velocity of Column 4
U4Init = np.array([[URexp[0],Umin,Umin,Umin,Umin,URexp[0],Umin,URexp[0],URexp[0],URexp[0],URexp[0],URexp[0]]])

for i in range(1,NData):
    U1Init = np.vstack((U1Init, np.array([[URexp[i],URaffexp[i],Umin,Umin,Umin,URexp[i],UExtexp[i],URexp[i],URexp[i],URexp[i],URexp[i],URexp[i]]])))
    U2Init = np.vstack((U2Init, np.array([[URexp[i],URaffexp[i],Umin,Umin,Umin,URexp[i],Umin,URexp[i],URexp[i],URexp[i],URexp[i],URexp[i]]]))) 
    U3Init = np.vstack((U3Init, np.array([[URexp[i],URaffexp[i],UFeedRaffexp[i],UFeedRaffexp[i],UFeedRaffexp[i],URexp[i],Umin,URexp[i],URexp[i],URexp[i],URexp[i],URexp[i]]])))
    U4Init = np.vstack((U4Init, np.array([[URexp[i],Umin,Umin,Umin,Umin,URexp[i],Umin,URexp[i],URexp[i],URexp[i],URexp[i],URexp[i]]])))

Colum = 4  #Number of columns
Nfet = 12  #Number of Finite Elements in Time 
Nfex = 30  #Number of Finite Elements in Space
ColL = 1.0 #Column length
UB = 8.0   #
LB = 0.1   #

Ncstr = 4  #Number of CSTR

# fill the initial parameter value
KapInit = [1.0, 1.0]
KInit = [1.0, 1.0]
bInit = [1.0, 1.0]

# -------------------------------------------------------------------
# Regularization Parameter
# -------------------------------------------------------------------

Rho0 = [4.0, 8.0]
Rho = {}
for i in range(2):
    Rho[i+1] = Rho0[i]

# fill the reliable parameter value
L2Kap = [1.0, 1.0]
L2K = [1.0, 1.0]
L2b = [1.0, 1.0]