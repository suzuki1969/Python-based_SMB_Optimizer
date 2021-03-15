# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kensuke Suzuki in 2020
# The default parameter setting in this file is the same as study by Sreedhar and Kawajiri
# Bibliographic information: Sreedhar B, Kawajiri Y. Multi-column chromatographic process development using simulated moving bed superstructure and simultaneous optimization - model correction framework. Chem Eng Sci. 2014;116:428-441. doi:10.1016/j.ces.2014.05.004

# -------------------------------------------------------------------
# Initial Values for Optimized Variables
# -------------------------------------------------------------------

U1Init = 1.460/3600 # Superfacial Velocity [m/s] in Colmun 1
U2Init = 1.053/3600 # Superfacial Velocity [m/s] in Colmun 2
U3Init = 1.256/3600 # Superfacial Velocity [m/s] in Colmun 3
U4Init = 0.951/3600 # Superfacial Velocity [m/s] in Colmun 4
StepTimeInit = 880.0 # Step Time

# -------------------------------------------------------------------
# Constraints Values for Control Variables
# -------------------------------------------------------------------

UB = 10.0 # Upper bound for Superfacial Velocity in Columns and Desorbent Velocity [m/hr]
LB = 0.1 # Lower bound for Throughput [m/hr]

ExtractRecMin0 = [0.0, 80] # Extract Recovery Bounds [A, B]
ExtractPurMin0 = [0.0, 90] # Extract Purity Bounds [A, B]
RaffinateRecMin0 = [0.0, 0.0] # Raffinate Recovery Bounds [A, B]
RaffinatePurMin0 = [0.0, 0.0] # Raffinate Purity Bounds [A, B]

# -------------------------------------------------------------------
# Discritization
# -------------------------------------------------------------------

Nfet = 5  # Number of Finite Elements in Time 
Nfex = 40  # Number of Finite Elements in Space per column
NCP = 3  # Number of Collocation Points
DScheme = 'CENTRAL' # 'BACKWARD' or 'CENTRAL'; Discretize scheme of space

# -------------------------------------------------------------------
# Details of subjected SMB
# -------------------------------------------------------------------

NComp = 2  # Number of components
NZone = 4 # This program only works only Four-Zone
Colum = 4  # Number of columns
ColL = 0.452  # Column length
eb = 0.375  # Void fraction
CFeed0 = [175, 175] # Feed concentration [A, B]

# -------------------------------------------------------------------
# LDF Model difinition 
# -------------------------------------------------------------------

Phase = 'Liquid' # 'Liquid', 'Solid' : Select LDF model based on liquid or solid phase 

Isotype = "Henry" # "Henry", "Langmuir" or "anti-Langmuir"

Axial_D = True # boolean : If True, axial dispersion implemented in LDF model

DV = True # Boolean : If True, Dead Volume is taken into account
Ncstr = 30  # Number of CSTRs

# -------------------------------------------------------------------
# Power Feed
# -------------------------------------------------------------------

PowerFeed = True # boolean : Power Feed
HT_Const = True # boolean : If True, PowerFeed works in the same time intervals 

# -------------------------------------------------------------------
# Model Parameter Values
# -------------------------------------------------------------------

# fill parameter values
Kap0 = [0.0047, 0.0083] # Overall Mass Transfer Coefficient [A, B]
H0 = [0.301, 0.531] # Henry's Constant [A, B]
b0 = [0.000634, 0.000240] # Affinity Coefficient [A, B]
Dax0 = [1.0e-8, 1.0e-8] # Axial Dispersion Coefficient [A, B]
DeadVolume = 0.05*Colum*ColL # Dead Volume

# -------------------------------------------------------------------
# Initializtion plots
# -------------------------------------------------------------------

ini_plots = False # boolean : Output internal concentration profiles at initilization

# -------------------------------------------------------------------
# Initial List Difinition
# -------------------------------------------------------------------

CFeed = {}
Kap = {}
H = {}
b = {}
Dax = {}
ExtractRecMin = {}
ExtractPurMin = {}
RaffinateRecMin = {}
RaffinatePurMin = {}

for i in range(NComp):
    CFeed[i+1] = CFeed0[i]
    Kap[i+1] = Kap0[i]
    H[i+1] = H0[i]
    b[i+1] = b0[i]
    Dax[i+1] = Dax0[i]
    ExtractRecMin[i+1] = ExtractRecMin0[i]
    ExtractPurMin[i+1] = ExtractPurMin0[i]
    RaffinateRecMin[i+1] = RaffinateRecMin0[i]
    RaffinatePurMin[i+1] = RaffinatePurMin0[i]
