# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kensuke Suzuki in 2020
# The default parameter setting in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

# -------------------------------------------------------------------
# Initial Values for Optimized Variables
# -------------------------------------------------------------------

U1Init = 0.002222222 # Superfacial Velocity [m/hr] in Colmun 1
U2Init = 0.001388889 # Superfacial Velocity [m/hr] in Colmun 2
U3Init = 0.001391667 # Superfacial Velocity [m/hr] in Colmun 3
U4Init = 0.000972222 # Superfacial Velocity [m/hr] in Colmun 4
StepTimeInit = 1140.0 # Step Time

# -------------------------------------------------------------------
# Constraints Values for Control Variables
# -------------------------------------------------------------------

UB = 8.0 # Upper bound for Superfacial Velocity in Columns and Desorbent Velocity
LB = 0.1 # Lower bound for Throughput

ExtractRecMin0 = [0.0, 80] # Extract Recovery Bounds [A, B]
ExtractPurMin0 = [0.0, 90] # Extract Purity Bounds [A, B]
RaffinateRecMin0 = [0.0, 0.0] # Raffinate Recovery Bounds [A, B]
RaffinatePurMin0 = [0.0, 0.0] # Raffinate Purity Bounds [A, B]

# -------------------------------------------------------------------
# Discritization
# -------------------------------------------------------------------

Nfet = 5  # Number of Finite Elements in Time 
Nfex = 40  # Number of Finite Elements in Space
NCP = 3  # Number of Collocation Points
DScheme = 'CENTRAL' # 'BACKWARD' or 'CENTRAL'; Discretize scheme of space

# -------------------------------------------------------------------
# Details of subjected ASMB
# -------------------------------------------------------------------

NComp = 2  # Number of components
NZone = 4 # This program only works only Four-Zone
Colum = 4  # Number of columns
ColL = 2.0  # Column length
eb = 0.389  # Void fraction
CFeed0 = [50.0, 50.0] # Feed concentration [A, B]
Isotype = "Henry" # "Henry", "Langmuir" or "anti-Langmuir"

# -------------------------------------------------------------------
# Dead volume
# -------------------------------------------------------------------

DV = 'yes' # 'yes', 'no' : Dead Volume
Ncstr = 4  # Number of CSTRs
DeadVolume = 0.05*Colum*ColL # Value of Dead Volume

# -------------------------------------------------------------------
# Power Feed
# -------------------------------------------------------------------

PowerFeed = 'yes' # 'yes', 'no' : Power Feed
HT_Const = 'no' # 'yes', 'no':If yes, PowerFeed works in the same time intervals 

# -------------------------------------------------------------------
# Model Parameter Values
# -------------------------------------------------------------------

# fill parameter values
Kap0 = [0.00684, 0.00684] # Overall Mass Transfer Coefficient [A, B]
K0 = [0.518, 0.743] # Henry's Constant [A, B]
b0 = [0.00002, 0.00002] # Affinity Coefficient [A, B]

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

for i in range(NComp):
    CFeed[i+1] = CFeed0[i]
    Kap[i+1] = Kap0[i]
    K[i+1] = K0[i]
    b[i+1] = b0[i]
    ExtractRecMin[i+1] = ExtractRecMin0[i]
    ExtractPurMin[i+1] = ExtractPurMin0[i]
    RaffinateRecMin[i+1] = RaffinateRecMin0[i]
    RaffinatePurMin[i+1] = RaffinatePurMin0[i]


# -------------------------------------------------------------------
# Show Details of SMB Model and Descritization
# -------------------------------------------------------------------

print(f"\n----- Discritization Scheme for Axcial Coordinate: {DScheme} -----")
print(f"\n----- Equilibrium Isotherm: {Isotype} -----")