# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# ASMB model was constructed by Hideki Harada and Kensuke Suzuki in 2020
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

import numpy as np

# -------------------------------------------------------------------
# Discritization
# -------------------------------------------------------------------

Nfet = 12  # Number of Finite Elements in Time 
Nfex = 30  # Number of Finite Elements in Space
NCP = 3  # Number of Collocation Points
DScheme = 'CENTRAL' # 'BACKWARD' or 'CENTRAL'; Discretize scheme of space

# -------------------------------------------------------------------
# Details of subjected ASMB
# -------------------------------------------------------------------

NComp = 2  # Number of components
NZone = 4 # This program only works only Four-Zone
Colum = 4  # Number of columns
ColL = 1.5  # Column length
eb = 0.4  # Void fraction
Isotype = "anti-Langmuir" # "Henry", "Langmuir" or "anti-Langmuir"

# -------------------------------------------------------------------
# Dead volume
# -------------------------------------------------------------------

Ncstr = 4  # Number of CSTRs

# -------------------------------------------------------------------
# Initial Values for Parameter Estimation
# -------------------------------------------------------------------

# fill the initial parameter value
KapInit = [1.0, 1.0]
KInit = [1.0, 1.0]
bInit = [1.0, 1.0]

# -------------------------------------------------------------------
# Tikhonov Regularization
# -------------------------------------------------------------------

Rho0 = [3.0, 8.0]
Rho = {}
for i in range(2):
    Rho[i+1] = Rho0[i]

# fill a reliable parameter value
L2Kap = [1.0, 1.0]
L2K = [1.0, 1.0]
L2b = [1.0, 1.0]

# -------------------------------------------------------------------
# Loading Data in CSV Files 
# -------------------------------------------------------------------

# fill the number of dataset in CSV files 
NData = 30 # Number Of Dataset

########################################

filename = "SubStepTime.csv"
ALLHT = np.loadtxt(filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

if NData == 1:
    HT1 = ALLHT[0]
    HT2 = ALLHT[1]
    HT3 = ALLHT[2]
    HT4 = ALLHT[3]
    HT5 = ALLHT[4]
    HT6 = ALLHT[5]
    StepTimeInit = [sum(ALLHT)]
else:
    HT1 = ALLHT[0,0:NData]
    HT2 = ALLHT[1,0:NData]
    HT3 = ALLHT[2,0:NData]
    HT4 = ALLHT[3,0:NData]
    HT5 = ALLHT[4,0:NData]
    HT6 = ALLHT[5,0:NData]
    StepTimeInit = [sum(ALLHT[0:6,i]) for i in range(NData)]

filename = "SubStepVelocity.csv"
ALLU = np.loadtxt(filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

if NData == 1:
    UICexp  = ALLU[0]
    URaffexp = ALLU[1]
    UFeedRaffexp = ALLU[2]
    UExtexp = ALLU[3]
else:
    UICexp  = ALLU[0,0:NData]
    URaffexp = ALLU[1,0:NData]
    UFeedRaffexp = ALLU[2,0:NData]
    UExtexp = ALLU[3,0:NData]
Umin = 0.000001

filename = "ProductConcentration.csv"
ALLConc = np.loadtxt(filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

CFeed = {}
CExtract = {}
CRaffinate = {}

if NData == 1:
    CFeed[1,1] = ALLConc[0]
    CFeed[1,2] = ALLConc[1]
    CRaffinate[1,1] = ALLConc[2]
    CRaffinate[1,2] = ALLConc[3]
    CExtract[1,1] = ALLConc[4]
    CExtract[1,2] = ALLConc[5]
else:
    Concentration_of_Feed_A = ALLConc[0,0:NData]
    Concentration_of_Feed_B = ALLConc[1,0:NData]
    Concentration_of_Raffinate_A = ALLConc[2,0:NData]
    Concentration_of_Raffinate_B = ALLConc[3,0:NData]
    Concentration_of_Extract_A = ALLConc[4,0:NData]
    Concentration_of_Extract_B = ALLConc[5,0:NData]

    for i in range(NData):
        CFeed[i+1,1] = Concentration_of_Feed_A[i]
        CExtract[i+1,1] = Concentration_of_Extract_A[i]
        CRaffinate[i+1,1] = Concentration_of_Raffinate_A[i]
        CFeed[i+1,2] = Concentration_of_Feed_B[i]
        CExtract[i+1,2] = Concentration_of_Extract_B[i]
        CRaffinate[i+1,2] = Concentration_of_Raffinate_B[i]

if NData == 1:
    HTRatio = Nfet*ALLHT/np.mean(StepTimeInit)
else:
    HTRatio = Nfet*np.mean(ALLHT, axis=1)/np.mean(StepTimeInit)
HTR_Int = [round(HTRatio[i]) for i in range(len(HTRatio))]

index = []
for i in range(len(HTRatio)):
    if HTR_Int[i] == 0:
        index.append(i)
        HTR_Int[i] = HTR_Int[i]+1
    elif  HTR_Int[i] == 1:
        index.append(i)
    else:
        pass

Val = sorted(enumerate([abs(HTRatio[i] - HTR_Int[i]) for i in range(len(HTRatio))]), reverse=True, key=lambda x:x[1])
for i in range(len(HTRatio)):
    if sum(HTR_Int) != Nfet:
        if sum(HTR_Int) < Nfet:
            HTR_Int[Val[i][0]] = HTR_Int[Val[i][0]] + 1
        else:
            if i in index:
                pass
            else:
                HTR_Int[Val[i][0]] = HTR_Int[Val[i][0]] - 1
    else:
        break
HTR = np.array(HTR_Int, dtype='int64')

# Time descritization of Feed-Raffinate sub-step should be finer
HTR[np.argmax(HTR)] = HTR[np.argmax(HTR)]-1
HTR[2] += 1

UInit = [[],[],[],[]]
for i in range(NData):
    if NData == 1:
         U = [[UICexp,URaffexp,Umin,UICexp,UExtexp,UICexp],
              [UICexp,URaffexp,Umin,UICexp,Umin,UICexp],
              [UICexp,URaffexp,UFeedRaffexp,UICexp,Umin,UICexp],
              [UICexp,Umin,Umin,UICexp,Umin,UICexp],]
    else:
        U = [[UICexp[i],URaffexp[i],Umin,UICexp[i],UExtexp[i],UICexp[i]],
             [UICexp[i],URaffexp[i],Umin,UICexp[i],Umin,UICexp[i]],
             [UICexp[i],URaffexp[i],UFeedRaffexp[i],UICexp[i],Umin,UICexp[i]],
             [UICexp[i],Umin,Umin,UICexp[i],Umin,UICexp[i]],]
    for j in range(Colum):
        UIte = []
        for k in range(len(HTRatio)):
            for l in range(HTR[k]):
                UIte.append(U[j][k])
        UInit[j].append(UIte)

U1Init = np.array(UInit[0][0]) # Initial values of fluid velocity in Zone 1
U2Init = np.array(UInit[1][0]) # Initial values of fluid velocity in Zone 2
U3Init = np.array(UInit[2][0]) # Initial values of fluid velocity in Zone 3
U4Init = np.array(UInit[3][0]) # Initial values of fluid velocity in Zone 4 

for i in range(1,NData):
    U1Init = np.vstack((U1Init, np.array(UInit[0][i])))
    U2Init = np.vstack((U2Init, np.array(UInit[1][i]))) 
    U3Init = np.vstack((U3Init, np.array(UInit[2][i])))
    U4Init = np.vstack((U4Init, np.array(UInit[3][i])))

HTI = NCP*HTR # Using for "HT_Fix_rule" and "Ucon"