# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kosuke Sekishita and Kensuke Suzuki in 2021
# The method used in this file is the same as study by Kawajiri and Biegler 
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.
# The superstructuer approach used in this script comes from Sreedhar and Kawajiri
# Bibliographic information: Sreedhar B, Kawajiri Y. Multi-column chromatographic process development using simulated moving bed superstructure and simultaneous optimization - model correction framework. Chem Eng Sci. 2014;116:428-441

import numpy as np

# -------------------------------------------------------------------
# Model Parameter Values
# -------------------------------------------------------------------

Kap0 = [0.00684, 0.00684] # Overall Mass Transfer Coefficient [A, B]
K0 = [0.518, 0.743] # Henry's Constant [A, B]
b0 = [0.0001, 0.0002] # Affinity Coefficient [A, B]

Isotype = "Henry"    # "Henry", "Langmuir" or "anti-Langmuir"

# -------------------------------------------------------------------
# Details of subjected SMB
# -------------------------------------------------------------------

NComp = 2  # Number of components
NZone = 4  # This program only works only Four-Zone
Colum = 4  # Number of columns

ColL = 2.0  # Column length [m]
ColD = 0.015  # Column diameter [m]
eb   = 0.389  # Void fraction [-]

CFeed0 = [175.0, 175.0] # Feed concentration [A, B] [g/L]

"""
ColL = 1.0    # Column length [m]
ColD = 0.01   # Column diameter [m]
eb = 0.3858   # Void fraction [-]

CFeed0 = [100.0, 100.0] # Feed concentration [A, B] [g/L]
"""

# -------------------------------------------------------------------
# Initial Values for Optimized Variables
# -------------------------------------------------------------------
U = np.empty((5,2))

"""
# Description as volumetric flow rates
Qrec = 0.1395 *60 #[mL/min]
QD   = 0.0414 *60 #[mL/min]
QE   = 0.0348 *60 #[mL/min]
QF   = 0.0200 *60 #[mL/min]

Q1 = Qrec
Q2 = Q1 - QE
Q3 = Q2 + QF
Q4 = Qrec - QD

U1Init = Q1 * 10**(-6) / (60 * np.pi*(1.3*10**(-2))**2) #[m/s]
U2Init = Q2 * 10**(-6) / (60 * np.pi*(1.3*10**(-2))**2) #[m/s]
U3Init = Q3 * 10**(-6) / (60 * np.pi*(1.3*10**(-2))**2) #[m/s]
U4Init = Q4 * 10**(-6) / (60 * np.pi*(1.3*10**(-2))**2) #[m/s]

UNInit = [U1Init, U2Init, U3Init, U4Init]
"""

# Description as flow velocity
Copy='''
U[1,1] = 7.600 [m/hr]
U[2,1] = 7.400 [m/hr]
U[3,1] = 7.500 [m/hr]
U[4,1] = 7.200 [m/hr]
'''.splitlines()[1:5]
UNInit = list(map(lambda x : float(x.split()[2])/3600, Copy))

StepTimeInit = 900.0 # Step Time [s]

# -------------------------------------------------------------------
# Constraints Values for Control Variables
# -------------------------------------------------------------------

#Q_UB = 3.0 # Upper bound for Flow Rate in Column [mL/min]
#UB = (Q_UB*(1e-6*60))/((ColD**2)*np.pi/4) # Upper bound for Superfacial Velocity in Columns and Desorbent Velocity [m/hr]
UB = 10.04 #[m/hr]

Q_LB = 0.0 #Lower bound for Feed Flow Rate [mL/min]
LB = (Q_LB*(1e-6*60))/((ColD**2)*np.pi/4) # Lower bound for Throughput [m/hr]

dS = 0.0 #the safety factor for purity [%]

ExtractRecMin0    = [0.0, 80.0]  # Extract Recovery Bounds [A, B]
ExtractPurMin0    = [0.0, 90.0] # Extract Purity Bounds [A, B]
RaffinateRecMin0  = [0.0, 0.0]  # Raffinate Recovery Bounds [A, B]
RaffinatePurMin0  = [0.0, 0.0] # Raffinate Purity Bounds [A, B]
ConcentrationMin0 = [0.0, 0.0]  # Concentration Bounds [A of Raf, B of Ext]

# -------------------------------------------------------------------
# Dead volume
# -------------------------------------------------------------------

DV = 'no' # 'yes', 'no' : Dead Volume
Ncstr = 20  # Number of CSTRs
DeadVolume = 0.0225*Colum*ColL # Value of Dead Volume

# -------------------------------------------------------------------
# Power Feed
# -------------------------------------------------------------------

PowerFeed = 'yes' # 'yes', 'no' : Power Feed
HT_Const  = 'no' # 'yes', 'no':If yes, PowerFeed works in the same time intervals

# -------------------------------------------------------------------
# Discritization
# -------------------------------------------------------------------

Nfet = 5  # Number of Finite Elements in Time = Number of Substep
Nfex = 40  # Number of Finite Elements in Space
NCP  = 3  # Number of Collocation Points

DScheme = 'CENTRAL' # 'BACKWARD' or 'CENTRAL'; Discretize scheme of space

# -------------------------------------------------------------------
# Directory Name Difinition
# -------------------------------------------------------------------

name = 'SMB'

if PowerFeed == 'yes':
    name += '_PowerFeed'
    print("\n----- Power Feed is implemented -----")
    if DV == 'yes':
        name += '_DeadVolume'
        print("\n----- Dead Volume is implemented -----")
    elif DV == 'no':
        pass
    else:
        print("ERROR: Confirm DeadVolume definition")
elif PowerFeed == 'no':
    HT_Const = 'no'
    if DV == 'yes':
        name += '_DeadVolume'
        print("\n----- Dead Volume is implemented -----")
    elif DV == 'no':
        pass
    else:
        print("ERROR: Confirm DeadVolume definition")
else:
    print("ERROR: Confirm PowerFeed definition")

# -------------------------------------------------------------------
# Initial List Difinition
# -------------------------------------------------------------------

CFeed = {}
Kap = {}
K = {}
b = {}
ExtractRecMin = {}
ExtractPurMin = {}
RaffinateRecMin = {}
RaffinatePurMin = {}
ConcentrationMin = {}

for i in range(NComp):
    CFeed[i+1] = CFeed0[i]
    Kap[i+1] = Kap0[i]
    K[i+1] = K0[i]
    b[i+1] = b0[i]
    ExtractRecMin[i+1] = ExtractRecMin0[i]
    ExtractPurMin[i+1] = ExtractPurMin0[i]
    RaffinateRecMin[i+1] = RaffinateRecMin0[i]
    RaffinatePurMin[i+1] = RaffinatePurMin0[i]
    ConcentrationMin[i+1] = ConcentrationMin0[i]

# -------------------------------------------------------------------
# Show Details of SMB Model and Descritization
# -------------------------------------------------------------------
Notification = [f"----- Discritization Scheme for Axcial Coordinate : {DScheme} -----", f"----- Equilibrium Isotherm : {Isotype} -----"]