# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kensuke Suzuki in 2020
# The method used in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

# ==================================================================
# Pyomo script file
# ==================================================================

from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import SolverFactory
from SMB_Model import m
from Initdata_SMB import *
from GaussRadauQuadrature import lglnodes
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import datetime
import sys
import os

dirname00 = f"Result_of_{name}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
os.makedirs(dirname00, exist_ok=True)

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open(f"{dirname00}Output_ParamEst_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.txt", "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        pass
    
sys.stdout = Logger()

sys.stderr = Logger()

Nfex = Nfex+1
Zone = int(Colum/NZone)

# -------------------------------------------------------------------
# Internal Velocities Constraint
# -------------------------------------------------------------------

def UMaxBound_rule(m,j,k):
	return m.U[j,k]*3600 <= UB

m.UMaxBound = Constraint(m.Col, m.t, rule = UMaxBound_rule)

# -------------------------------------------------------------------
# Set Constraint List
# -------------------------------------------------------------------

m.U_constraints=ConstraintList()
m.HT_constraints=ConstraintList()

# -------------------------------------------------------------------
# Create Instance
# -------------------------------------------------------------------

instance = m.create_instance()

# -------------------------------------------------------------------
# Discretize time and space using collocation and finite difference, respectively
# -------------------------------------------------------------------
#Discretize model using Radau Collocation method 
discretizer = TransformationFactory('dae.collocation')
#Discretize model using finite difference method 
discretizer2 = TransformationFactory('dae.finite_difference')

#nfe=number of finite elements
#ncp=number of collocation points within each finite element, (scheme=OO, default is ‘LAGRANGE-RADAU’)
discretizer.apply_to(instance,nfe=Nfet, ncp=value(instance.NCP) , wrt=instance.t, scheme='LAGRANGE-RADAU')

discretizer2.apply_to(instance,nfe=Nfex-1, wrt=instance.x, scheme='CENTRAL')

if DV == 'yes':
    discretizer2.apply_to(instance,nfe=Ncstr, wrt=instance.xRe, scheme='BACKWARD')
elif DV == 'no':
    pass
else:
    print("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# Difinition of  Objective Function
# -------------------------------------------------------------------

def obj_expr(m):
    return 3600*(sum(m.HT[k]*m.UF[j,k]/m.NCP for j in m.Col for k in m.t if k != 0))

instance.obj = Objective(rule = obj_expr, sense = maximize)

# -------------------------------------------------------------------
# Throughput Constraint
# -------------------------------------------------------------------

def MinThroughput_constraint(m):
	return 3600*sum(m.HT[k]*m.UF[j,k]/m.NCP for j in m.Col for k in m.t if k != 0) >= LB

instance.MinThroughput = Constraint(rule = MinThroughput_constraint)

# -------------------------------------------------------------------
# Desorbent Constraint
# -------------------------------------------------------------------

def UDMaxBound_constraint(m):
	return 3600*sum(m.HT[k]*m.UD[j,k]/m.NCP for j in m.Col for k in m.t if k != 0) <= UB

instance.UDMaxBound = Constraint(rule = UDMaxBound_constraint)

# -------------------------------------------------------------------
# Fixing variables using Build Action
# -------------------------------------------------------------------

# ------------------------------------------------------
# Fixing internal fluid velocities in each column
# ------------------------------------------------------

# Fluid velocities in column are the same at all NCP in a sub-step

def U1Fix_rule(m):
    for k in m.t:
        m.U[1,k].fix(U1Init)

instance.U1Fix = BuildAction(rule = U1Fix_rule)

def U2Fix_rule(m):
    for k in m.t:
        m.U[Zone+1,k].fix(U2Init)

instance.U2Fix = BuildAction(rule = U2Fix_rule)

def U3Fix_rule(m):
    for k in m.t:
        m.U[2*Zone+1,k].fix(U3Init)

instance.U3Fix = BuildAction(rule = U3Fix_rule)

def U4Fix_rule(m):
    for k in m.t:
        m.U[3*Zone+1,k].fix(U4Init)

instance.U4Fix = BuildAction(rule = U4Fix_rule)

# ------------------------------------------------------
# Fixing velocities from unused inlet/putlet ports at 0
# ------------------------------------------------------

# Feed is only supplied from inlet of first colmun in Zone 3
def UF_Fix_rule(m):
    for i in m.Col:
        if i!=(2*Zone+1):
            for k in m.t:
                m.UF[i,k].fix(0)

instance.UF_Fix = BuildAction(rule = UF_Fix_rule)

# Extract is only collected from outlet of last colmun in Zone 1
def UE_Fix_rule(m):
    for i in m.Col:
        if i!=(1*Zone):
            for k in m.t:
                m.UE[i,k].fix(0)

instance.UE_Fix = BuildAction(rule = UE_Fix_rule)

# Desorbent is only supplied from inlet of first colmun in Zone 1
def UD_Fix_rule(m):
    for i in m.Col:
        if i!=(0*Zone+1):
            for k in m.t:
                m.UD[i,k].fix(0)

instance.UD_Fix = BuildAction(rule = UD_Fix_rule)

# Raffinate is only collected from outlet of last colmun in Zone 3
def UR_Fix_rule(m):
    for i in m.Col:
        if i!=(3*Zone):
            for k in m.t:
                m.UR[i,k].fix(0)

instance.UR_Fix = BuildAction(rule = UR_Fix_rule)

def StepTimeFix_rule(m):
	m.StepTime.fix(StepTimeInit)

instance.StepTimeFix = BuildAction(rule = StepTimeFix_rule)

# ------------------------------------------------------
# Fixing HT (Finite elements along time)
# ------------------------------------------------------

def HT_Fix_rule(m):
    for i in m.t:
        m.HT[i].fix(1/Nfet)

instance.HT_Fix = BuildAction(rule = HT_Fix_rule)

def HTsum_rule(m):
    return (sum(m.HT[i] for i in m.t if i != 0) == m.NCP)

instance.HTSum = Constraint(rule = HTsum_rule)

for i in range(0,Nfet):
    for j in range(2,value(instance.NCP)+1):
        instance.HT_constraints.add(instance.HT[instance.t[i*instance.NCP+j]] == instance.HT[instance.t[i*instance.NCP+j+1]])
instance.HT_constraints.add(instance.HT[instance.t[1]] == instance.HT[instance.t[2]])

# -------------------------------------------------------------------
# Calculating weights of Gauss-Radau quadrature
# -------------------------------------------------------------------

w=lglnodes(value(instance.NCP))

def Ac_init(m,cp):
    return w[cp-1]/2.0
instance.Ac = Param(instance.CP,initialize=Ac_init,default=0)

#---------------------------------------------------------------------
# Additional Constraint
#---------------------------------------------------------------------

def _intCE(m, i,j):
    return sum( sum(m.C[i,j,m.t[k*m.NCP+1+l],m.L]*m.Ac[l]for l in m.CP)*m.UE[j,m.t[k*m.NCP+2]]*m.HT[m.t[k*m.NCP+2]] for k in range(0,Nfet)) == m.intCE[i,j]
instance.intCE0 = Constraint(instance.Comp, instance.Col, rule=_intCE)

def _intCR(m, i,j):
    return sum( sum(m.C[i,j,m.t[k*m.NCP+1+l],m.L]*m.Ac[l] for l in m.CP)*m.UR[j,m.t[k*m.NCP+2]]*m.HT[m.t[k*m.NCP+2]] for k in range(0,Nfet)) == m.intCR[i,j]
instance.intCR0 = Constraint(instance.Comp, instance.Col, rule=_intCR)

# -------------------------------------------------------------------
# Activate if CENTRAL scheme is selected
# -------------------------------------------------------------------

# Describing last axial derivate as a backward finite-difference method with 2nd-order 
if DScheme == 'CENTRAL':
    def C_AxialDerivativeConstraintEnd_rule(m,i,j,k):
        return m.dCdx[i,j,k,value(m.L)]==(3*m.C[i,j,k,value(m.L)]-4*m.C[i,j,k,m.x[Nfex-1]]+m.C[i,j,k,m.x[Nfex-2]])/(2*m.L/(Nfex-1)) 
    instance.C_AxialDerivativeConstraintEnd=Constraint(instance.Comp, instance.Col, instance.t,  rule=C_AxialDerivativeConstraintEnd_rule)
else:
    pass

# -------------------------------------------------------------------
# Activate if CENTRAL scheme is selected
# -------------------------------------------------------------------

def ExtractRecoveryConstraint_rule(m):
	return sum(m.intCE[2,j] for j in m.Col) >= (m.ExtractRecMin[2]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[2]/m.NCP for j in m.Col for k in m.t if k != 0)

def RaffinateRecoveryConstraint_rule(m):
	return sum(m.intCR[1,j] for j in m.Col) >= (m.RaffinateRecMin[1]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[1]/m.NCP for j in m.Col for k in m.t if k != 0)

instance.ExtractRecoveryConstraint = Constraint(rule=ExtractRecoveryConstraint_rule)
instance.RaffinateRecoveryConstraint = Constraint(rule=RaffinateRecoveryConstraint_rule)

# -------------------------------------------------------------------
# Over All Mass Balance Constraint
# -------------------------------------------------------------------

def OverallBalance_rule(m,k):
	return sum(m.UR[j,k] for j in m.Col ) == sum(m.UD[j,k]+m.UF[j,k]-m.UE[j,k] for j in m.Col )

instance.OverallBalance=Constraint(instance.t,rule=OverallBalance_rule)

# -------------------------------------------------------------------
# Activate and Deactivate Constraints
# -------------------------------------------------------------------

instance.ExtractRecoveryConstraint.deactivate()
instance.ExtractPurityConstraint.deactivate()
instance.RaffinateRecoveryConstraint.deactivate()
instance.RaffinatePurityConstraint.deactivate()
instance.UMaxBound.deactivate()
instance.MinThroughput.deactivate()
instance.UDMaxBound.deactivate()
instance.HTSum.deactivate()
instance.U_constraints.deactivate()
instance.HT_constraints.deactivate()

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------

opt = SolverFactory('ipopt')
opt.options['mu_init'] = 1e-3
opt.options['ma57_pivtol'] = 1e-8
opt.options['max_iter'] = 5000
opt.options['linear_system_scaling'] = 'mc19'

# -------------------------------------------------------------------
# Solve successively
# -------------------------------------------------------------------

instance.Lvar.fix()
instance.preprocess()
results = opt.solve(instance, tee=True)

# file = open('tmp.txt','w')
# instance.display(ostream=file)
# file.close()

# -------------------------------------------------------------------
# Display Results
# -------------------------------------------------------------------

print("\n-------------------Initial Condition Solution----------------------\n")

for i in instance.Comp:
    for j in instance.Col:
        print(value(instance.intCE[i,j]))

for i in instance.Comp:
    for j in instance.Col:
        print(value(instance.intCR[i,j]))

for i in instance.Comp:
	print("Comp A product recovery %s : \t %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

for i in instance.Comp:
	print("Comp B product purity %s : \t %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.intCE[n,j] for n in instance.Comp for j in instance.Col))))

for i in instance.Comp:
	print("Comp A product recovery %s : \t %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col )/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

for i in instance.Comp:
	print("Comp B product purity %s : \t %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col)/sum(instance.intCR[n,j] for n in instance.Comp for j in instance.Col ))))

print('Control Varialbes: \n')
if PowerFeed == 'yes':
    print('U')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("U[%s,%s] = %s [m/hr]" %(i,k,value(3600*instance.U[i,k])))

    print('UD')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UD[%s] = %s [m/hr]" %(i,value(3600*instance.UD[i,k])))

    print('UE')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UE[%s] = %s [m/hr]" %(i,value(3600*instance.UE[i,k])))

    print('UF')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UF[%s] = %s [m/hr]" %(i,value(3600*instance.UF[i,k])))
    print('UR')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UR[%s] = %s [m/hr]" %(i,value(3600*instance.UR[i,k])))

elif PowerFeed == 'no':
    print('U')
    for i in sorted(instance.Col):
        print ("U[%s,1] = %s [m/hr]" %(i,value(3600*instance.U[i,1])))

    print('UD')
    for i in sorted(instance.Col):
        print ("UD[%s] = %s [m/hr]" %(i,value(3600*instance.UD[i,1])))

    print('UE')
    for i in sorted(instance.Col):
        print ("UE[%s] = %s [m/hr]" %(i,value(3600*instance.UE[i,1])))

    print('UF')
    for i in sorted(instance.Col):
        print ("UF[%s] = %s [m/hr]" %(i,value(3600*instance.UF[i,1])))

    print('UR')
    for i in sorted(instance.Col):
        print ("UR[%s] = %s [m/hr]" %(i,value(3600*instance.UR[i,1])))

else:
    print("ERROR: Confirm PowerFeed definition")

print("Step Time %s [s]" %(value(instance.StepTime)))

print("Comp A conc in extract = %s" %(value((sum(instance.intCE[1,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Comp B conc in extract = %s" %(value((sum(instance.intCE[2,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Comp A conc in raffinate = %s" %(value((sum(instance.intCR[1,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Comp B conc in raffinate = %s" %(value((sum(instance.intCR[2,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))


# -------------------------------------------------------------------
# Dropping constraints 
# -------------------------------------------------------------------

for k in instance.t:
    instance.U[1,k].free()
    instance.U[Zone+1,k].free()
    instance.U[2*Zone+1,k].free()
    instance.U[3*Zone+1,k].free()
    instance.HT[k].free()

instance.StepTime.free()

# -------------------------------------------------------------------
# Setting flow rate of ASMB
# -------------------------------------------------------------------

if PowerFeed == 'yes':
    pass
elif PowerFeed == 'no':
    for i in range(1,Nfet):
        instance.HT_constraints.add(instance.HT[instance.t[2]] == instance.HT[instance.t[i*instance.NCP+2]])
    
    for i in range(1,Nfet):
        instance.U_constraints.add(instance.U[Zone+1,instance.t[2]] == instance.U[Zone+1,instance.t[i*instance.NCP+2]])

    for i in range(1,Nfet):
        instance.U_constraints.add(instance.U[1,instance.t[2]] == instance.U[1,instance.t[i*instance.NCP+2]])

    for i in range(1,Nfet):
        instance.U_constraints.add(instance.U[2*Zone+1,instance.t[2]] == instance.U[2*Zone+1,instance.t[i*instance.NCP+2]])

    for i in range(1,Nfet):
        instance.U_constraints.add(instance.U[3*Zone+1,instance.t[2]] == instance.U[3*Zone+1,instance.t[i*instance.NCP+2]])
else:
    print("ERROR: Confirm PowerFeed definition")

if HT_Const == 'yes':
    for i in range(1,Nfet):
        instance.HT_constraints.add(instance.HT[instance.t[2]] == instance.HT[instance.t[i*instance.NCP+2]])
elif HT_Const == 'no':
    pass
else:
    print("ERROR: Confirm HT_COnst definition")

for i in range(0,Nfet):
    for j in range(2,value(instance.NCP)+1):
        instance.U_constraints.add(instance.U[Zone+1,instance.t[i*instance.NCP+j]] == instance.U[Zone+1,instance.t[i*instance.NCP+j+1]])
instance.U_constraints.add(instance.U[Zone+1,instance.t[1]] == instance.U[Zone+1,instance.t[2]])

for i in range(0,Nfet):
    for j in range(2,value(instance.NCP)+1):
        instance.U_constraints.add(instance.U[1,instance.t[i*instance.NCP+j]] == instance.U[1,instance.t[i*instance.NCP+j+1]])
instance.U_constraints.add(instance.U[1,instance.t[1]] == instance.U[1,instance.t[2]])

for i in range(0,Nfet):
    for j in range(2,value(instance.NCP)+1):
        instance.U_constraints.add(instance.U[2*Zone+1,instance.t[i*instance.NCP+j]] == instance.U[2*Zone+1,instance.t[i*instance.NCP+j+1]])
instance.U_constraints.add(instance.U[2*Zone+1,instance.t[1]] == instance.U[2*Zone+1,instance.t[2]])

for i in range(0,Nfet):
    for j in range(2,value(instance.NCP)+1):
        instance.U_constraints.add(instance.U[3*Zone+1,instance.t[i*instance.NCP+j]] == instance.U[3*Zone+1,instance.t[i*instance.NCP+j+1]])
instance.U_constraints.add(instance.U[3*Zone+1,instance.t[1]] == instance.U[3*Zone+1,instance.t[2]])

# -------------------------------------------------------------------
# Activating constraints
# -------------------------------------------------------------------

instance.ExtractRecoveryConstraint.activate()
instance.ExtractPurityConstraint.activate()
instance.RaffinateRecoveryConstraint.activate()
instance.RaffinatePurityConstraint.activate()
instance.UMaxBound.activate()
instance.MinThroughput.activate()
instance.UDMaxBound.activate()

instance.HTSum.activate()
instance.U_constraints.activate()
instance.HT_constraints.activate()

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------
opt.options['mu_init'] = 1e-3
opt.options['ma57_pivtol'] = 1e-8
opt.options['max_iter'] = 5000
opt.options['linear_system_scaling'] = 'mc19'
opt.options['linear_solver']='ma27'

instance.preprocess()
results = opt.solve(instance, tee=True)

print("\n-------------------Optimzed Solution----------------------\n")

if PowerFeed == 'yes':
    print("\n----- Power Feed is implemented -----")
    if DV == 'yes':
        print("\n----- Dead Volume is implemented -----")
    elif DV == 'no':
        pass
    else:
        print("ERROR: Confirm DeadVolume definition")
elif PowerFeed == 'no':
    if DV == 'yes':
        print("\n----- Dead Volume is implemented -----")
    elif DV == 'no':
        pass
    else:
        print("ERROR: Confirm DeadVolume definition")
else:
    print("ERROR: Confirm PowerFeed definition")

print(f"\n----- Discritization Scheme for Axcial Coordinate :{DScheme} -----")
print(f"\n----- Equilibrium Isotherm :{Isotype} -----")

for i in instance.Comp:
	print("\nComp A product recovery %s : \t %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

for i in instance.Comp:
	print("Comp B product purity %s : \t %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.intCE[n,j] for n in instance.Comp for j in instance.Col))))

for i in instance.Comp:
	print("Comp A product recovery %s : \t %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col )/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

for i in instance.Comp:
	print("Comp B product purity %s : \t %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col)/sum(instance.intCR[n,j] for n in instance.Comp for j in instance.Col ))))

print('Control Varialbes: \n')
if PowerFeed == 'yes':
    print('U')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("U[%s,%s] = %s [m/hr]" %(i,k,value(3600*instance.U[i,k])))

    print('UD')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UD[%s] = %s [m/hr]" %(i,value(3600*instance.UD[i,k])))

    print('UE')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UE[%s] = %s [m/hr]" %(i,value(3600*instance.UE[i,k])))

    print('UF')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UF[%s] = %s [m/hr]" %(i,value(3600*instance.UF[i,k])))
    print('UR')
    for i in sorted(instance.Col):
        for k in instance.t:
            print ("UR[%s] = %s [m/hr]" %(i,value(3600*instance.UR[i,k])))
    t = []   
    for i in sorted(instance.t):
        t.append(value(instance.StepTime)*i)
    print(f"\nTime span in Steptime [s]: {t}")

elif PowerFeed == 'no':
    print('U')
    for i in sorted(instance.Col):
        print ("U[%s,1] = %s [m/hr]" %(i,value(3600*instance.U[i,1])))

    print('UD')
    for i in sorted(instance.Col):
        print ("UD[%s] = %s [m/hr]" %(i,value(3600*instance.UD[i,1])))

    print('UE')
    for i in sorted(instance.Col):
        print ("UE[%s] = %s [m/hr]" %(i,value(3600*instance.UE[i,1])))

    print('UF')
    for i in sorted(instance.Col):
        print ("UF[%s] = %s [m/hr]" %(i,value(3600*instance.UF[i,1])))

    print('UR')
    for i in sorted(instance.Col):
        print ("UR[%s] = %s [m/hr]" %(i,value(3600*instance.UR[i,1])))

else:
    print("ERROR: Confirm PowerFeed definition")

print("\nStep Time %s [s]" %(value(instance.StepTime)))

print("\nThroughput: %s [m/hr]" %(value(3600*sum(instance.UF[j,k]*instance.HT[k]/instance.NCP for k in instance.t if k != 0 for j in instance.Col))))

print("\nDesorbent: %s [m/hr]" %(value(3600*sum(instance.UD[j,k]*instance.HT[k]/instance.NCP for k in instance.t if k != 0 for j in instance.Col))))

print("\nD/F ratio: %s" %(value((sum(instance.UD[j,k]*instance.HT[k]/instance.NCP for k in instance.t for j in instance.Col))/(sum(instance.UF[j,k]*instance.HT[k]/instance.NCP for k in instance.t if k != 0 for j in instance.Col)))))

print("Comp A conc in extract = %s" %(value((sum(instance.intCE[1,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Comp B conc in extract = %s" %(value((sum(instance.intCE[2,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Comp A conc in raffinate = %s" %(value((sum(instance.intCR[1,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Comp B conc in raffinate = %s" %(value((sum(instance.intCR[2,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Optimization Completed")

# -------------------------------------------------------------------
# Plot graphs
# -------------------------------------------------------------------

dirname01 = dirname00 + f"Plotfile_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
os.makedirs(dirname01, exist_ok=True)

X = []
CA = []
CB = []

for i in sorted(instance.Col):
    for j in sorted(instance.x):
       CA.append((value(instance.C[1, i, 1.0,j]))/value(instance.CF[1]))
       CB.append((value(instance.C[2, i, 1.0,j]))/value(instance.CF[2]))
       X.append(j+(i-1)*ColL)
plt.figure()
plt.plot(X, CA, 'o-', label='Comp A')
plt.plot(X, CB, 'x-', label='Comp B')
plt.legend()
plt.xlim(0, ColL*Colum)
plt.ylim(-1.0e-2, 1.5)
plt.xlabel('x [m]')
plt.ylabel('C/C(feed) [-]')
plt.title(f"Concentration profile")
filename = dirname01 + f"concentration_profile"
plt.savefig(filename, dpi=500)
plt.close()

dirname02 = dirname00 + f"Animationfile_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
os.makedirs(dirname02, exist_ok=True)

C1all=[]
C2all=[]

for k in sorted(instance.t):
    C1 = []
    C2 = []
    for i in sorted(instance.Col):
        for j in sorted(instance.x):
            C1.append((value(instance.C[1, i, k ,j]))/value(instance.CF[1]))
            C2.append((value(instance.C[2, i, k ,j]))/value(instance.CF[2]))
    C1all.append(C1)
    C2all.append(C2)

C1all = np.array(C1all)
C2all = np.array(C2all)

C1_ani = C1all
C2_ani = C2all

C1all = np.hstack([C1all[:,3*Nfex:],C1all[:,:3*Nfex]])
C2all = np.hstack([C2all[:,3*Nfex:],C2all[:,:3*Nfex]])
C1_ani = np.vstack([C1_ani,C1all])
C2_ani = np.vstack([C2_ani,C2all])

C1all = np.hstack([C1all[:,3*Nfex:],C1all[:,:3*Nfex]])
C2all = np.hstack([C2all[:,3*Nfex:],C2all[:,:3*Nfex]])
C1_ani = np.vstack([C1_ani,C1all])
C2_ani = np.vstack([C2_ani,C2all])

C1all = np.hstack([C1all[:,3*Nfex:],C1all[:,:3*Nfex]])
C2all = np.hstack([C2all[:,3*Nfex:],C2all[:,:3*Nfex]])
C1_ani = np.vstack([C1_ani,C1all])
C2_ani = np.vstack([C2_ani,C2all])

fig = plt.figure(figsize = (8, 6))

def update(i,u1,u2,x,fig_title):
    plt.cla()
    plt.ylim(-1.0e-2, 1.5)
    plt.xlim(0, ColL*Colum)
    plt.plot(x,u1[i,:], 'o-', label='Comp A')
    plt.plot(x,u2[i,:], 'x-', label='Comp B')
    plt.xlabel('x [m]')
    plt.ylabel('C/C(feed) [-]')
    plt.title(f"SMB Concentration profile")
    plt.legend()
    plt.close()

nframe = Nfet*value(instance.NCP)+1
nframe = nframe*4
ani = animation.FuncAnimation(fig,update,fargs = (C1_ani,C2_ani,X,'Wave motion'),interval=100,frames=nframe)
filename = dirname02 + f"SMB_Concentration_Animation.gif"
ani.save(filename, writer="pillow", dpi=500)

if PowerFeed == 'yes':
    dirname03 = dirname00 + f"PowerFeed_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
    os.makedirs(dirname03, exist_ok=True)
    
    dirname04 = dirname00 + f"Internal_Profile_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
    os.makedirs(dirname04, exist_ok=True)

    t = []
    Colum
    U_Col = [[] for i in range(Colum)]
    U_four = [[] for i in range(NZone)]
    U_name = ["Desorbent","Extract","Feed","Raffinate"]
    U_color = ["C2","C3","C4","C5"]

    for i in sorted(instance.t):
        U_four[0].append(value(3600*instance.UD[0*Zone+1,i]))
        U_four[1].append(value(3600*instance.UE[1*Zone,i]))
        U_four[2].append(value(3600*instance.UF[2*Zone+1,i]))
        U_four[3].append(value(3600*instance.UR[3*Zone,i]))
        t.append(value(instance.StepTime)*i)

    for i in range(NZone):
        plt.figure()
        plt.plot(t, U_four[i], 'x-', color=U_color[i], label=U_name[i])
        plt.legend()
        plt.xlim(0, value(instance.StepTime))
        plt.ylim(-1.0e-2, 1.2*max(U_four[i]))
        plt.xlabel('Step Time [s]')
        plt.ylabel(f"{U_name[i]} [m/hr]")
        plt.title(f"Profile of {U_name[i]}")
        filename = dirname03 + f"{U_name[i]}_profile"
        plt.savefig(filename, dpi=500)
        plt.close()

    for j in sorted(instance.Col):
        for i in sorted(instance.t):
            U_Col[j-1].append(value(3600*instance.U[j,i]))

    for i in range(Colum):
        plt.figure()
        plt.plot(t, U_Col[i], 'x-', color='k')
        plt.xlim(0, value(instance.StepTime))
        plt.ylim(-1.0e-2, 1.2*max(U_Col[i]))
        plt.xlabel('Step Time [s]')
        plt.ylabel(f"Internal Superfacial Velocity [m/hr]")
        plt.title(f"Superfacial Velocity in Column {i+1}")
        filename = dirname04 + f"Column_{i+1}_profile"
        plt.savefig(filename, dpi=500)
        plt.close()
elif PowerFeed == 'no':
    pass
else:
    print("ERROR: Confirm PowerFeed definition")
