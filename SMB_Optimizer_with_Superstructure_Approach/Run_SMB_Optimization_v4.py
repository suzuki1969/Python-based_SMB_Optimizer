# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kosuke Sekishita and Kensuke Suzuki in 2021
# The method used in this file is the same as study by Kawajiri and Biegler 
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.
# The superstructuer approach used in this script comes from Sreedhar and Kawajiri
# Bibliographic information: Sreedhar B, Kawajiri Y. Multi-column chromatographic process development using simulated moving bed superstructure and simultaneous optimization - model correction framework. Chem Eng Sci. 2014;116:428-441

# ==================================================================
# Pyomo script file
# ==================================================================

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math
import datetime
import sys
import os

from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import SolverFactory

from GaussRadauQuadrature import lglnodes
from SMB_Model_v1 import m
from Initdata_SMB_v1 import *

Con = []
NameList = []

# -------------------------------------------------------------------
# File Specification
# -------------------------------------------------------------------

dirname00 = f"{name}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
os.makedirs(dirname00, exist_ok=True)

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open(f"{dirname00}Output_Optimization_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.txt", "a")
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
# Set Constraint List
# -------------------------------------------------------------------

m.UF_constraint=ConstraintList()
m.UD_constraint=ConstraintList()
m.UE_constraint=ConstraintList()
m.UR_constraint=ConstraintList()
m.HT_constraint=ConstraintList()

# -------------------------------------------------------------------
# Slack variable
# -------------------------------------------------------------------

m.Slack = Var(initialize=1e-2, within=PositiveReals)

# -------------------------------------------------------------------
# Create Instance
# -------------------------------------------------------------------

instance = m.create_instance()

# -------------------------------------------------------------------
# Discretize time and space using collocation and finite difference, respectively
# -------------------------------------------------------------------

discretizer = TransformationFactory('dae.collocation')        #Discretize model using Radau Collocation method
discretizer2 = TransformationFactory('dae.finite_difference') #Discretize model using finite difference method

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

from Objective_Function.ThMax_with_Slack import *
exec(compile(ObjectiveFunction, "", "exec"))
Notification += Notif

# -------------------------------------------------------------------
# Add Constraint Files
# -------------------------------------------------------------------

from Constraint.FeConst_ver0 import *
from Constraint.DeConst_slack_ver0 import *
from Constraint.ThConst_ver0 import *
from Constraint.RaConst_ver0 import *
from Constraint.ExConst_ver0 import *

Con += [FeConst, DeConst, ThConst, RaConst, ExConst]

for i in Con:
    exec(compile(i, "", "exec"))
    NameList += Name
    Notification += Notif
print()
print("\n\n".join(Notification),"\n")

# ------------------------------------------------------
# Fixing internal fluid velocities in each column
# ------------------------------------------------------
@instance.BuildAction()
def UallFix(m):
    [m.U[i*Zone+1, k].fix(UNInit[i]) for i in range(4) for k in m.t]

# ------------------------------------------------------
# Fixing velocities from unused inlet/outlet ports at 0
# ------------------------------------------------------

# Desorbent is only supplied from inlet of first colmun in Zone 1
@instance.BuildAction()
def UD_Fix(m):
    [m.UD[i,k].fix(0) for i in m.Col if i!=0*Zone+1 for k in m.t]

# Raffinate is only collected from outlet of last colmun in Zone 3
@instance.BuildAction()
def UR_Fix(m):
    [m.UR[i,k].fix(0) for i in m.Col if i!=3*Zone for k in m.t]

# Feed is only supplied from inlet of first colmun in Zone 3
@instance.BuildAction()
def UF_Fix(m):
    [m.UF[i,k].fix(0) for i in m.Col if i!=2*Zone+1 for k in m.t]

# Extract is only collected from outlet of last colmun in Zone 1
@instance.BuildAction()
def UE_Fix(m):
    [m.UE[i,k].fix(0) for i in m.Col if i!=1*Zone for k in m.t]

@instance.BuildAction()
def StepTimeFix(m):
	m.StepTime.fix(StepTimeInit)

# ------------------------------------------------------
# Fixing HT (Finite elements along time)
# ------------------------------------------------------
@instance.BuildAction()
def HT_Fix(m):
    [m.HT[i].fix(1/Nfet) for i in m.t]

@instance.Constraint()
def HTSum(m):
    return (sum(m.HT[i] for i in m.t if i != 0) == m.NCP)

[instance.HT_constraint.add(instance.HT[instance.t[i*instance.NCP+j]] == instance.HT[instance.t[i*instance.NCP+j+1]]) \
 for i in range(0,Nfet) for j in range(2,value(instance.NCP)+1)]
instance.HT_constraint.add(instance.HT[instance.t[1]] == instance.HT[instance.t[2]])

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

@instance.Constraint(instance.Comp, instance.Col)
def intCE0(m,i,j):
    return sum( sum(m.C[i,j,m.t[k*m.NCP+1+l],m.L]*m.Ac[l] for l in m.CP)*m.UE[j,m.t[k*m.NCP+2]]*m.HT[m.t[k*m.NCP+2]] for k in range(0,Nfet)) == m.intCE[i,j]
@instance.Constraint(instance.Comp, instance.Col)
def intCR0(m,i,j):
    return sum( sum(m.C[i,j,m.t[k*m.NCP+1+l],m.L]*m.Ac[l] for l in m.CP)*m.UR[j,m.t[k*m.NCP+2]]*m.HT[m.t[k*m.NCP+2]] for k in range(0,Nfet)) == m.intCR[i,j]

# -------------------------------------------------------------------
# Activate if CENTRAL scheme is selected
# -------------------------------------------------------------------

# Describing last axial derivate as a backward finite-difference method with 2nd-order
if DScheme == 'CENTRAL':
    @instance.Constraint(instance.Comp, instance.Col, instance.t)
    def C_AxialDerivativeConstraintEnd(m,i,j,k):
        return m.dCdx[i,j,k,value(m.L)]==(3*m.C[i,j,k,value(m.L)]-4*m.C[i,j,k,m.x[Nfex-1]]+m.C[i,j,k,m.x[Nfex-2]])/(2*m.L/(Nfex-1))
else:
    pass

# -------------------------------------------------------------------
# Activate if CENTRAL scheme is selected
# -------------------------------------------------------------------

@instance.Constraint()
def ExtractRecoveryConstraint(m):
	return sum(m.intCE[2,j] for j in m.Col) >= (m.ExtractRecMin[2]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[2]/m.NCP for j in m.Col for k in m.t if k != 0)
@instance.Constraint()
def RaffinateRecoveryConstraint(m):
	return sum(m.intCR[1,j] for j in m.Col) >= (m.RaffinateRecMin[1]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[1]/m.NCP for j in m.Col for k in m.t if k != 0)

# -------------------------------------------------------------------
# Over All Mass Balance Constraint
# -------------------------------------------------------------------

@instance.Constraint(instance.t)
def OverallBalance(m,k):
	return sum(m.UR[j,k] for j in m.Col ) == sum(m.UD[j,k]+m.UF[j,k]-m.UE[j,k] for j in m.Col )

# -------------------------------------------------------------------
# Deactivate Constraints
# -------------------------------------------------------------------

instance.ExtractRecoveryConstraint.deactivate()
instance.ExtractPurityConstraint.deactivate()
instance.RaffinateRecoveryConstraint.deactivate()
instance.RaffinatePurityConstraint.deactivate()
instance.HTSum.deactivate()
instance.HT_constraint.deactivate()
[eval(f"instance.{i}.deactivate()") for i in NameList]

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------

opt = SolverFactory('ipopt')
opt.options['mu_init'] = 1e-3
opt.options['ma57_pivtol'] = 1e-8 #収束性
opt.options['max_iter'] = 5000
# opt.options['linear_system_scaling'] = 'mc19'
opt.options['linear_solver'] = 'ma97'

# -------------------------------------------------------------------
# Solve successively
# -------------------------------------------------------------------

instance.Lvar.fix()
instance.preprocess()
opt.solve(instance, tee=True)

# file = open('tmp.txt','w')
# instance.display(ostream=file)
# file.close()

# -------------------------------------------------------------------
# Display Results
# -------------------------------------------------------------------

print("\n-------------------Initial Condition Solution----------------------\n")

for i in instance.Comp:
    print(f"【Concentration {i} in Exract】")
    for j in instance.Col:
        print(f"Colmn{j} : ",value(instance.intCE[i,j]))

for i in instance.Comp:
    print(f"【Concentration {i} in Raffinate】")
    for j in instance.Col:
        print(f"Colmn{j} : ",value(instance.intCR[i,j]))

print("\n【Recovary】")
for i in instance.Comp:
	print("%s in Extract   : %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))
for i in instance.Comp:
	print("%s in Raffinate : %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col )/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

print("\n【Purity】")
for i in instance.Comp:
	print("%s in Extract   : %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.intCE[n,j] for n in instance.Comp for j in instance.Col))))
for i in instance.Comp:
	print("%s in Raffinate : %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col)/sum(instance.intCR[n,j] for n in instance.Comp for j in instance.Col ))))

print('\n【Control Varialbes】')
if PowerFeed == 'yes':
    print('U [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" U[%s,%s] = %s" %(i,k,value(3600*instance.U[i,k])))

    print('UD [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UD[%s]  = %s" %(i,value(3600*instance.UD[i,k])))

    print('UE [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UE[%s]  = %s" %(i,value(3600*instance.UE[i,k])))

    print('UF [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UF[%s]  = %s" %(i,value(3600*instance.UF[i,k])))
    print('UR [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UR[%s]  = %s" %(i,value(3600*instance.UR[i,k])))

elif PowerFeed == 'no':
    print('U [m/hr]')
    for i in sorted(instance.Col):
        print (" U[%s,1] = %s" %(i,value(3600*instance.U[i,1])))

    print('UD [m/hr]')
    for i in sorted(instance.Col):
        print (" UD[%s]  = %s" %(i,value(3600*instance.UD[i,1])))

    print('UE [m/hr]')
    for i in sorted(instance.Col):
        print (" UE[%s]  = %s" %(i,value(3600*instance.UE[i,1])))

    print('UF [m/hr]')
    for i in sorted(instance.Col):
        print (" UF[%s]  = %s" %(i,value(3600*instance.UF[i,1])))

    print('UR [m/hr]')
    for i in sorted(instance.Col):
        print (" UR[%s]  = %s" %(i,value(3600*instance.UR[i,1])))

else:
    print("ERROR: Confirm PowerFeed definition")

print("\nStep Time  : %s [s]" %(value(instance.StepTime)))

print("\n【Concentration】")

print("A conc in extract   = %s [g/L]" %(value((sum(instance.intCE[1,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("B conc in extract   = %s [g/L]" %(value((sum(instance.intCE[2,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("A conc in raffinate = %s [g/L]" %(value((sum(instance.intCR[1,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("B conc in raffinate = %s [g/L]\n" %(value((sum(instance.intCR[2,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

# -------------------------------------------------------------------
# Dropping constraints
# -------------------------------------------------------------------

[instance.U[i*Zone+1,k].free() for i in range(4) for k in instance.t]
[instance.HT[k].free() for k in instance.t]
instance.StepTime.free()

[instance.UE[i,k].free() for i in instance.Col if i!=1*Zone for k in instance.t]
[instance.UR[i,k].free() for i in instance.Col if i!=3*Zone for k in instance.t]
[instance.UF[i,k].free() for i in instance.Col if i!=2*Zone+1 for k in instance.t]

# -------------------------------------------------------------------
# Setting flow rate of ASMB
# -------------------------------------------------------------------

if PowerFeed == 'yes':
    pass
elif PowerFeed == 'no':
    [instance.HT_constraint.add(instance.HT[instance.t[2]] == instance.HT[instance.t[i*instance.NCP+2]]) for i in range(1,Nfet)]

    [instance.UE_constraint.add(instance.U[Zone+1,instance.t[2]] == instance.U[Zone+1,instance.t[i*instance.NCP+2]]) for i in range(1,Nfet)]

    [instance.UD_constraint.add(instance.U[1,instance.t[2]] == instance.U[1,instance.t[i*instance.NCP+2]]) for i in range(1,Nfet)]

    [instance.UF_constraint.add(instance.U[2*Zone+1,instance.t[2]] == instance.U[2*Zone+1,instance.t[i*instance.NCP+2]]) for i in range(1,Nfet)]

    [instance.UR_constraint.add(instance.U[3*Zone+1,instance.t[2]] == instance.U[3*Zone+1,instance.t[i*instance.NCP+2]]) for i in range(1,Nfet)]

    [instance.UF_constraint.add(instance.UF[h,instance.t[2]] == instance.UF[h,instance.t[i*instance.NCP+2]]) for h in range(1, 4*Zone+1) for i in range(1,Nfet)]

    [instance.UE_constraint.add(instance.UE[h,instance.t[2]] == instance.UE[h,instance.t[i*instance.NCP+2]]) for h in range(1, 4*Zone+1) for i in range(1,Nfet)]

    [instance.UR_constraint.add(instance.UR[4*Zone,instance.t[2]] == instance.UR[4*Zone,instance.t[i*instance.NCP+2]]) for i in range(1,Nfet)]

else:
    print("ERROR: Confirm PowerFeed definition")

if HT_Const == 'yes':
    [instance.HT_constraint.add(instance.HT[instance.t[2]] == instance.HT[instance.t[i*instance.NCP+2]]) for i in range(1,Nfet)]
elif HT_Const == 'no':
    pass
else:
    print("ERROR: Confirm HT_COnst definition")

# for i in range(0,Nfet):
#     for j in range(2,value(instance.NCP)+1):
#         instance.UE_constraint.add(instance.U[Zone+1,instance.t[i*instance.NCP+j]] == instance.U[Zone+1,instance.t[i*instance.NCP+j+1]])

[instance.UE_constraint.add(instance.U[Zone+1,instance.t[i*instance.NCP+j]] == instance.U[Zone+1,instance.t[i*instance.NCP+j+1]]) for i in range(0,Nfet) for j in range(2,value(instance.NCP)+1)]
instance.UE_constraint.add(instance.U[Zone+1,instance.t[1]] == instance.U[Zone+1,instance.t[2]])

[instance.UD_constraint.add(instance.U[1,instance.t[i*instance.NCP+j]] == instance.U[1,instance.t[i*instance.NCP+j+1]]) for i in range(0,Nfet) for j in range(2,value(instance.NCP)+1)]
instance.UD_constraint.add(instance.U[1,instance.t[1]] == instance.U[1,instance.t[2]])

[instance.UF_constraint.add(instance.U[2*Zone+1,instance.t[i*instance.NCP+j]] == instance.U[2*Zone+1,instance.t[i*instance.NCP+j+1]]) for i in range(0,Nfet) for j in range(2,value(instance.NCP)+1)]
instance.UF_constraint.add(instance.U[2*Zone+1,instance.t[1]] == instance.U[2*Zone+1,instance.t[2]])

[instance.UR_constraint.add(instance.U[3*Zone+1,instance.t[i*instance.NCP+j]] == instance.U[3*Zone+1,instance.t[i*instance.NCP+j+1]]) for i in range(0,Nfet) for j in range(2,value(instance.NCP)+1)]
instance.UR_constraint.add(instance.U[3*Zone+1,instance.t[1]] == instance.U[3*Zone+1,instance.t[2]])

for h in range(1, 4*Zone+1):
    for i in range(0,Nfet):
        for j in range(2,value(instance.NCP)+1):
            instance.UF_constraint.add(instance.UF[h,instance.t[i*instance.NCP+j]] == instance.UF[h,instance.t[i*instance.NCP+j+1]])
    instance.UF_constraint.add(instance.UF[h,instance.t[1]] == instance.UF[h,instance.t[2]])

for h in range(1, 4*Zone+1):
    for i in range(0,Nfet):
        for j in range(2,value(instance.NCP)+1):
            instance.UE_constraint.add(instance.UE[h,instance.t[i*instance.NCP+j]] == instance.UE[h,instance.t[i*instance.NCP+j+1]])
    instance.UE_constraint.add(instance.UE[h,instance.t[1]] == instance.UE[h,instance.t[2]])

[instance.UR_constraint.add(instance.UR[4*Zone,instance.t[i*instance.NCP+j]] == instance.UR[4*Zone,instance.t[i*instance.NCP+j+1]]) for i in range(0,Nfet) for j in range(2,value(instance.NCP)+1)]
instance.UR_constraint.add(instance.UR[4*Zone,instance.t[1]] == instance.UR[4*Zone,instance.t[2]])

# -------------------------------------------------------------------
# Activating constraints
# -------------------------------------------------------------------

instance.ExtractRecoveryConstraint.activate()
instance.ExtractPurityConstraint.activate()
instance.RaffinateRecoveryConstraint.activate()
instance.RaffinatePurityConstraint.activate()
instance.HTSum.activate()
instance.HT_constraint.activate()
instance.UE_constraint.activate()
instance.UD_constraint.activate()
instance.UF_constraint.activate()
instance.UR_constraint.activate()
[eval(f"instance.{i}.activate()") for i in NameList]

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------
opt.options['mu_init'] = 1e-3
# opt.options['ma57_pivtol'] = 1e-8
opt.options['max_iter'] = 5000
# opt.options['linear_system_scaling'] = 'mc19'
opt.options['linear_solver'] = 'ma97'
# opt.options['tol'] = 1e-5

instance.preprocess()
opt.solve(instance, tee=True)

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

print(f"\n----- Discritization Scheme for Axcial Coordinate : {DScheme} -----")
print(f"\n----- Equilibrium Isotherm : {Isotype} -----\n")
print("\n\n".join(Notification))

print("\n【Recovary】")
for i in instance.Comp:
	print("%s in Extract   : %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))
for i in instance.Comp:
	print("%s in Raffinate : %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col )/sum(instance.HT[k]*instance.UF[j,k]*instance.CF[i]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

print("\n【Purity】")
for i in instance.Comp:
	print("%s in Extract   : %s" %(i, value(100*sum(instance.intCE[i,j] for j in instance.Col)/sum(instance.intCE[n,j] for n in instance.Comp for j in instance.Col))))
for i in instance.Comp:
	print("%s in Raffinate : %s" %(i, value(100*sum(instance.intCR[i,j] for j in instance.Col)/sum(instance.intCR[n,j] for n in instance.Comp for j in instance.Col ))))

print('\n【Control Varialbes】')
if PowerFeed == 'yes':
    print('U [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" U[%s,%s] = %s" %(i,k,value(3600*instance.U[i,k])))

    print('UD [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UD[%s]  = %s" %(i,value(3600*instance.UD[i,k])))

    print('UE [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UE[%s]  = %s" %(i,value(3600*instance.UE[i,k])))

    print('UF [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UF[%s]  = %s" %(i,value(3600*instance.UF[i,k])))
    print('UR [m/hr]')
    for i in sorted(instance.Col):
        for k in instance.t:
            print (" UR[%s]  = %s" %(i,value(3600*instance.UR[i,k])))
    t = []
    for i in sorted(instance.t):
        t.append(value(instance.StepTime)*i)
    print(f"\nTime span in Steptime [s] : {t}\n")

elif PowerFeed == 'no':
    print('U [m/hr]')
    for i in sorted(instance.Col):
        print (" U[%s,1] = %s" %(i,value(3600*instance.U[i,1])))

    print('UD [m/hr]')
    for i in sorted(instance.Col):
        print (" UD[%s]  = %s" %(i,value(3600*instance.UD[i,1])))

    print('UE [m/hr]')
    for i in sorted(instance.Col):
        print (" UE[%s]  = %s" %(i,value(3600*instance.UE[i,1])))

    print('UF [m/hr]')
    for i in sorted(instance.Col):
        print (" UF[%s]  = %s" %(i,value(3600*instance.UF[i,1])))

    print('UR [m/hr]')
    for i in sorted(instance.Col):
        print (" UR[%s]  = %s" %(i,value(3600*instance.UR[i,1])))

else:
    print("ERROR: Confirm PowerFeed definition")

print("\nStep Time  : %s [s]" %(value(instance.StepTime)))

print("\nThroughput : %s [m/hr]" %(value(3600*sum(instance.UF[j,k]*instance.HT[k]/instance.NCP for k in instance.t if k != 0 for j in instance.Col))))

print("\nDesorbent  : %s [m/hr]" %(value(3600*sum(instance.UD[j,k]*instance.HT[k]/instance.NCP for k in instance.t if k != 0 for j in instance.Col))))

print("\nD/F ratio  : %s [-]" %(value((sum(instance.UD[j,k]*instance.HT[k]/instance.NCP for k in instance.t for j in instance.Col))/(sum(instance.UF[j,k]*instance.HT[k]/instance.NCP for k in instance.t if k != 0 for j in instance.Col)))))

print("\n【Concentration】")

print("A in extract   = %s [g/L]" %(value((sum(instance.intCE[1,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("B in extract   = %s [g/L]" %(value((sum(instance.intCE[2,j] for j in instance.Col))/(sum(instance.HT[k]*instance.UE[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("A in raffinate = %s [g/L]" %(value((sum(instance.intCR[1,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("B in raffinate = %s [g/L]\n" %(value((sum(instance.intCR[2,j] for j in instance.Col ))/(sum(instance.HT[k]*instance.UR[j,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

print("Optimization Completed\n")

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
    #Colum
    U_Col = [[] for i in range(Colum)]
    U_four = [[] for i in range(NZone)]
    U_name = ["Desorbent","Extract","Feed","Raffinate"]
    U_color = ['gray',"blue","purple","red"]

#sekishita
    t = [value(instance.StepTime)*i for i in sorted(instance.t)]
    U_matrix = [[[]]*Colum for i in range(4)]
    for i in range(Colum):
        U_matrix[0][i] = list(map(lambda n: 3600*value(n), instance.UD[i+1,:]))
        U_matrix[1][i] = list(map(lambda n: 3600*value(n), instance.UE[i+1,:]))
        U_matrix[2][i] = list(map(lambda n: 3600*value(n), instance.UF[i+1,:]))
        U_matrix[3][i] = list(map(lambda n: 3600*value(n), instance.UR[i+1,:]))

    fig = plt.figure()
    plt.xlabel('Step Time [s]', labelpad=20)
    plt.xticks(ticks=[], labels=[])
    plt.ylabel("Superfacial Vel. of Nodes [m/hr]", labelpad=30)
    plt.yticks(ticks=[], labels=[])
    for i in range(Colum):
        ad = int(410+(i+1))
        ax = fig.add_subplot(ad)
        for j in range(4):
            ax.plot(t, U_matrix[j][i], 'x-', color=U_color[j], label=U_name[j])
        ax.set_xlim(0, value(instance.StepTime))
        ax.set_ylim(-1.0e-2, 1.2*max(list(map(max, U_matrix[:][i]))))
        ax.set_xlabel=('')
        if i != Colum-1:
            ax.tick_params(labelbottom=False)
        ax.tick_params(direction='in')
        ax.set_title(f"Column {i+1}", size=8, loc='right', pad=1)
    filename = dirname03 + f"Sup_Vel_of_Nodes"
    plt.legend(fontsize=5, loc='center right', bbox_to_anchor=(1.14,1.1))
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