# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kensuke Suzuki in 2020
# The default parameter setting in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.
# EVM formulation          : Zavala V, Biegler LT, Large-scale parameter estimation in low-density polyethylene tubular reactors. Industrial and Engineering Chemistry Research. 2006;45,23;7867-7881
# Parameter reference      : Sreedhar B, Kawajiri Y, Multi-column chromatographic process development using simulated moving bed superstructure and simultaneous optimization - Model correction framework, Chemical Engineering Science, 2014;116,428-441
# Tikhonov regularization  : Tie et al., Experimental evaluation of simulated moving bed reactor for transesterification reaction synthesis of glycol ether ester. Adsorption. 2019;25,4;795-807

import numpy as np
import sys

# -------------------------------------------------------------------
# Data
# -------------------------------------------------------------------

# fill the number of dataset in CSV files 
NData = 2 # Number Of Dataset

StpeTime_filename = "StepTimes.csv"
SubStepFlow_filename = "FlowRates.csv"
ProductConc_filename = "ProductConcentrations.csv"

# -------------------------------------------------------------------
# Discritization
# -------------------------------------------------------------------

Nfet = 5  # Number of Finite Elements in Time 
Nfex = 40  # Number of Finite Elements in Space
NCP = 3  # Number of Collocation Points
DScheme = 'CENTRAL' # 'BACKWARD' or 'CENTRAL'; Discretize scheme of space

# -------------------------------------------------------------------
# Details of subjected SMB
# -------------------------------------------------------------------

NComp = 2  # Number of components
NZone = 4 # This program only works only Four-Zone
Colum = 4  # Number of columns
ColL = 0.452  # Column length
eb = 0.375 # Void fraction

# -------------------------------------------------------------------
# LDF Model difinition 
# -------------------------------------------------------------------

Phase = 'Liquid' # 'Liquid', 'Solid' : Select LDF model based on liquid or solid phase 

Isotype = 'Henry' # "Henry", "Langmuir" or "anti-Langmuir"

Axial_D = False # boolean : If True, axial dispersion implemented in LDF model

DV = False # Boolean : If True, Dead Volume is taken into account
Ncstr = 4  # Number of CSTRs

# -------------------------------------------------------------------
# Power Feed
# -------------------------------------------------------------------

PowerFeed = False # boolean : Power Feed
HT_Const = True # boolean :If yes, PowerFeed works in the same time intervals 

# -------------------------------------------------------------------
# Model Parameter Values
# -------------------------------------------------------------------

# fill initial parameter values
Kapinit = [0.0047, 0.0083] # Overall Mass Transfer Coefficient [A, B]
Hinit = [0.301, 0.531] # Henry's Constant [A, B]
binit = [0.0, 0.0] # Affinity Coefficient [A, B]
Daxinit = [1.0e-5, 1.0e-5] # Axial Dispersion Coefficient [A, B]
LReinit = 0.05*ColL # Dead Volume Length

# -------------------------------------------------------------------
# Tikhonov Regularization and EVM estimation formulation
# -------------------------------------------------------------------

Tikhonov = True # boolean : Tikhonov Regularization
EVM = False # boolean : EVM Estimation Formulation
Regularization_Parameter = [1.0e1, 1.0e1] # [Tikhonov regularization parameter, EVM regularization parameter]
Nreg = [1,2]

# fill a reliable parameter value for Tikhonov Regularization
L2Kap = [0.0047, 0.0083] # Overall Mass Transfer Coefficient [A, B]
L2H = [0.301, 0.531] # Henry's Constant [A, B]
L2b = [0.0, 0.0] # Affinity Coefficient [A, B]
L2Dax = [1.0e-5, 1.0e-5] # Axial Dispersion Coefficient [A, B]
L2LRe = 0.05*Colum*ColL # Dead Volume Length

# -------------------------------------------------------------------
# Initializtion plots
# -------------------------------------------------------------------

ini_plots = True # boolean : Output internal concentration profiles at initilization

# -------------------------------------------------------------------
# For initializtion
# -------------------------------------------------------------------

ite_e = 3 #iteration for initilization (integer)

# -------------------------------------------------------------------
# Regularization and EVM estimation setting
# -------------------------------------------------------------------

if Tikhonov == True:
    if EVM == True:
        pass
    elif EVM == False:
        Regularization_Parameter[1] = 0
        pass
    else:
        sys.exit("ERROR: Confirm EVM definition")
elif Tikhonov == False:
    pass
else:
    sys.exit("ERROR: Confirm Tikhonov definition")

Rho = {}
for i in range(2):
    Rho[i+1] = Regularization_Parameter[i]

# -------------------------------------------------------------------
# Loading Data in CSV Files 
# -------------------------------------------------------------------

StepTimeInit = np.loadtxt(StpeTime_filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

if NData == 1:
    StepTimeInit = [StepTimeInit, 0]
else:
    pass

ALLU = np.loadtxt(SubStepFlow_filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

if NData == 1:
    UDinit  = [ALLU[0],0]
    UFinit = [ALLU[1],0]
    UEinit = [ALLU[2],0]
    UReinit = [ALLU[3],0]
else:
    UDinit  = ALLU[0,0:NData]
    UFinit = ALLU[1,0:NData]
    UEinit = ALLU[2,0:NData]
    UReinit = ALLU[3,0:NData]

ALLConc = np.loadtxt(ProductConc_filename, delimiter=",", skiprows=1, usecols=tuple(np.arange(1, NData+1)))

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
