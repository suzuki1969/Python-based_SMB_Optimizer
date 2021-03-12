# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# ASMB model was constructed by Hideki Harada and Kensuke Suzuki in 2020
# The method used in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

# ==================================================================
# Pyomo script file
# ==================================================================
from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import SolverFactory
from ParamEst_ASMB import m, data
from ParamEst_Initdata_ASMB import *
from GaussRadauQuadrature import lglnodes
import datetime
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import os
import numpy as np

dirname00 = f"Result_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
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

# Considering one step of SMB in this program

Nfex = Nfex+1
Zone = int(Colum/NZone)

m.Ucon = ConstraintList() # list of constraints on fluid velocity (used in parameter estimation step)

m.intCR = Var(m.Data, m.Comp, m.Col) # integrated concentration of Raffinate
m.intCE = Var(m.Data, m.Comp, m.Col) # integrated concentration of Extract

m.UofIC = Var(m.Data) # fluid velocity of Internal-Circulation-process in ASMB (used in parameter estimation step)
m.UofRaff = Var(m.Data) # fluid velocity of Raffinate-collected-process in ASMB
m.UofFeedRaff = Var(m.Data) # fluid velocity of Feed-injected-and-Raffinate-collected-process in ASMB
m.UofExt = Var(m.Data) # fluid velocity of Extract-collected-process in ASMB
m.Usmall = Param(initialize = 0.000001) # please reffer to the document. 

instance = m.create_instance() # For-loop works only after instance 

# -------------------------------------------------------------------
# Discretize time and space using collocation and finite difference, respectively
# -------------------------------------------------------------------
#Discretize model using Radau Collocation method 
discretizer = TransformationFactory('dae.collocation')
#Discretize model using finite difference method 
discretizer2 = TransformationFactory('dae.finite_difference')

#nfe=number of finite elements
#ncp=number of collocation points within each finite element, (scheme=OO, default is ‘LAGRANGE-RADAU’)
discretizer.apply_to(instance,nfe=Nfet,ncp=value(instance.NCP), wrt=instance.t)

discretizer2.apply_to(instance,nfe=Nfex-1, wrt=instance.x, scheme=DScheme)

discretizer2.apply_to(instance,nfe=Ncstr, wrt=instance.xRe, scheme='BACKWARD')

# -------------------------------------------------------------------
# Constraints of fluid velocity for each collocation point in finite elements
# -------------------------------------------------------------------

# In a sub-step, fluid velocity is the same to the next collocation point: "U[t[1]] == U[t[2]], U[t[2]] == U[t[3]], ... "

for d in sorted(instance.Data):
    for i in range(1,Colum+1):
        instance.Ucon.add(instance.U[d,(i-1)*Zone+1,instance.t[1]] == instance.U[d,(i-1)*Zone+1,instance.t[2]])
        for j in range(1,Nfet+1):
            for k in range(2,value(instance.NCP)+1):
                instance.Ucon.add(instance.U[d,(i-1)*Zone+1,instance.t[(j-1)*instance.NCP+k]] == instance.U[d,(i-1)*Zone+1,instance.t[(j-1)*instance.NCP+k+1]])

# -------------------------------------------------------------------
# Fixing variables using Build Action
# -------------------------------------------------------------------

# ------------------------------------------------------
# Fixing model parameter at initial value
# ------------------------------------------------------

def Kap_Fix_rule(m):
  for i in sorted(m.Comp):
      m.Kap[i].fix(KapInit[i-1])
instance.Kap_Fix = BuildAction(rule = Kap_Fix_rule)

def K_Fix_rule(m):
  for i in sorted(m.Comp):
      m.K[i].fix(KInit[i-1])
instance.K_Fix = BuildAction(rule = K_Fix_rule)

if Isotype == "Henry":
    pass
elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
    def b_Fix_rule(m):
        for i in sorted(m.Comp):
            m.b[i].fix(bInit[i-1])
    instance.b_Fix = BuildAction(rule = b_Fix_rule)
else:
    print("ERROR: Confirm Isotype definition")

# ------------------------------------------------------
# Fixing internal fluid velocities in each column
# ------------------------------------------------------

# Fluid velocities in column are the same at all NCP in a sub-step

def U1Fix_rule(m):
    for d in sorted(m.Data):
        for x in range(0*Zone+1, 1*Zone+1):
            m.U[d,x,0].fix(U1Init[d-1][0])
            for i in range(1,m.NCP+1):
                for j in range(0,Nfet):
                    m.U[d,x,m.t[m.NCP*j+i+1]].fix(U1Init[d-1][j])
instance.U1Fix = BuildAction(rule = U1Fix_rule)

def U2Fix_rule(m):
    for d in sorted(m.Data):
        for x in range(1*Zone+1, 2*Zone+1):
            m.U[d,x,0].fix(U2Init[d-1][0])
            for i in range(1,m.NCP+1):
                for j in range(0,Nfet):
                    m.U[d,x,m.t[m.NCP*j+i+1]].fix(U2Init[d-1][j])
instance.U2Fix = BuildAction(rule = U2Fix_rule)

def U3Fix_rule(m):
        for d in sorted(m.Data):
            for x in range(2*Zone+1, 3*Zone+1):
                m.U[d,x,0].fix(U3Init[d-1][0])
                for i in range(1,m.NCP+1):
                    for j in range(0,Nfet):
                        m.U[d,x,m.t[m.NCP*j+i+1]].fix(U3Init[d-1][j])
instance.U3Fix = BuildAction(rule = U3Fix_rule)

def U4Fix_rule(m):
    for d in sorted(m.Data):
        for x in range(3*Zone+1, 4*Zone+1):
            m.U[d,x,0].fix(U4Init[d-1][0])
            for i in range(1,m.NCP+1):
                for j in range(0,Nfet):
                    m.U[d,x,m.t[m.NCP*j+i+1]].fix(U4Init[d-1][j])
instance.U4Fix = BuildAction(rule = U4Fix_rule)

# ------------------------------------------------------
# Fixing velocities from unused inlet/putlet ports at 0
# ------------------------------------------------------

# Feed is only supplied from inlet of first colmun in Zone 3
def UF_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i!=(2*Zone+1):
                    m.UF[d,i,j].fix(0)
instance.UF_Fix = BuildAction(rule = UF_Fix_rule)

# Extract is only collected from outlet of last colmun in Zone 1
def UE_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i!=(1*Zone):
                    m.UE[d,i,j].fix(0)
instance.UE_Fix = BuildAction(rule = UE_Fix_rule)

# Desorbent is only supplied from inlet of first colmun in Zone 1
def UD_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i!=(0*Zone+1):
                    m.UD[d,i,j].fix(0)
instance.UD_Fix = BuildAction(rule = UD_Fix_rule)

# Raffinate is only collected from outlet of last colmun in Zone 3
def UR_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i!=(3*Zone):
                    m.UR[d,i,j].fix(0)
instance.UR_Fix = BuildAction(rule = UR_Fix_rule)

def StepTimeFix_rule(m):
    for d in sorted(m.Data):
        m.StepTime[d].fix(StepTimeInit[d-1])
instance.StepTimeFix = BuildAction(rule = StepTimeFix_rule)

# ------------------------------------------------------
# Fixing HT (Finite elements along time)
# ------------------------------------------------------

def HT_Fix_rule(m):
    for d in sorted(m.Data):
        for i in m.t:
            if i < value(m.t[1+sum(HTI[0:1])+1]):
                m.HT[d,i].fix(1.0*HT1[d-1]/(HTR[0]*StepTimeInit[d-1]))
            if value(m.t[1+sum(HTI[0:1])]) < i and i < value(m.t[1+sum(HTI[0:2])+1]):
                m.HT[d,i].fix(1.0*HT2[d-1]/(HTR[1]*StepTimeInit[d-1]))
            if value(m.t[1+sum(HTI[0:2])]) < i and i < value(m.t[1+sum(HTI[0:3])+1]):
                m.HT[d,i].fix(1.0*HT3[d-1]/(HTR[2]*StepTimeInit[d-1]))
            if value(m.t[1+sum(HTI[0:3])]) < i and i < value(m.t[1+sum(HTI[0:4])+1]):
                m.HT[d,i].fix(1.0*HT4[d-1]/(HTR[3]*StepTimeInit[d-1]))
            if value(m.t[1+sum(HTI[0:4])]) < i and i < value(m.t[1+sum(HTI[0:5])+1]):
                m.HT[d,i].fix(1.0*HT5[d-1]/(HTR[4]*StepTimeInit[d-1]))
            if value(m.t[1+sum(HTI[0:5])]) < i:
                m.HT[d,i].fix(1.0*HT6[d-1]/(HTR[5]*StepTimeInit[d-1]))
instance.HT_Fix = BuildAction(rule = HT_Fix_rule)

# -------------------------------------------------------------------
# Calculating weights of Gauss-Radau quadrature
# -------------------------------------------------------------------

w=lglnodes(value(instance.NCP))

def Ac_init(m,cp):
    return w[cp-1]/2.0
instance.Ac = Param(instance.CP,initialize=Ac_init,default=0)

# -------------------------------------------------------------------
# Defining calculation of Integrals with Gauss-Radau quadrature
# -------------------------------------------------------------------

def _intCR(m, d,i,j):
    return sum(m.UR[d,j,m.t[(k-1)*m.NCP+1+l]]*m.HT[d,m.t[(k-1)*m.NCP+1+l]]*(m.Ac[l]*m.C[d,i,j,m.t[1+(k-1)*m.NCP+l],value(m.L)]) for l in range(1,m.NCP+1) for k in range(1,Nfet+1)) == m.intCR[d,i,j]
instance.intCRrule = Constraint(instance.Data, instance.Comp, instance.Col,rule=_intCR)

def _intCE(m, d,i,j):
    return sum(m.UE[d,j,m.t[(k-1)*m.NCP+1+l]]*m.HT[d,m.t[(k-1)*m.NCP+1+l]]*(m.Ac[l]*m.C[d,i,j,m.t[1+(k-1)*m.NCP+l],value(m.L)]) for l in range(1,m.NCP+1) for k in range(1,Nfet+1)) == m.intCE[d,i,j]
instance.intCErule = Constraint(instance.Data, instance.Comp, instance.Col,rule=_intCE)

# -------------------------------------------------------------------
# Activate if CENTRAL scheme is selected
# -------------------------------------------------------------------

# Describing last axial derivate as a backward finite-difference method with 2nd-order 
if DScheme == 'CENTRAL':
    def C_AxialDerivativeConstraintEnd_rule(m,d,i,j,k):
        return m.dCdx[d,i,j,k,value(m.L)]==(3*m.C[d,i,j,k,value(m.L)]-4*m.C[d,i,j,k,m.x[Nfex-1]]+m.C[d,i,j,k,m.x[Nfex-2]])/(2*m.L/(Nfex-1)) 
    instance.C_AxialDerivativeConstraintEnd=Constraint(instance.Data, instance.Comp, instance.Col, instance.t,  rule=C_AxialDerivativeConstraintEnd_rule)
else:
    pass

# -------------------------------------------------------------------
# Objective Function
# -------------------------------------------------------------------

# Objective function 1 to maximize productivity, used for initial simulation in this program
def obj1_expr(m):
    return 3600*sum(sum(m.UF[1,j,k]*m.HT[1,k]/m.NCP for k in m.t if k != 0) for j in m.Col)
instance.obj1 = Objective(rule = obj1_expr, sense = maximize)

# Objective function 2 to estimate model parameters together with fluid velocities
def obj2_expr(m):
    if Isotype == "Henry":
        return (1/len(m.Data))*sum((m.CbarE[d,i] - m.averageCE[d,i])**2 + (m.CbarR[d,i] - m.averageCR[d,i])**2 for d in m.Data for i in m.Comp) + ((10**m.rho[2])/len(m.Data))*sum(((m.UofIC[d] - U1Init[d-1][0])/U1Init[d-1][0])**2 + ((m.UofRaff[d] - U1Init[d-1][1])/U1Init[d-1][1])**2 + ((m.UofFeedRaff[d] - U3Init[d-1][3])/U3Init[d-1][3])**2 + ((m.UofExt[d] - U1Init[d-1][6])/U1Init[d-1][6])**2 for d in m.Data) + (10**m.rho[1])*(sum(((L2K[i-1] - m.K[i])/L2K[i-1])**2 + ((L2Kap[i-1] - m.Kap[i])/L2Kap[i-1])**2 for i in m.Comp) + ((m.LRe - 0.05*Colum*ColL)/(0.05*Colum*ColL))**2)
    elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
        return (1/len(m.Data))*sum((m.CbarE[d,i] - m.averageCE[d,i])**2 + (m.CbarR[d,i] - m.averageCR[d,i])**2 for d in m.Data for i in m.Comp) + ((10**m.rho[2])/len(m.Data))*sum(((m.UofIC[d] - U1Init[d-1][0])/U1Init[d-1][0])**2 + ((m.UofRaff[d] - U1Init[d-1][1])/U1Init[d-1][1])**2 + ((m.UofFeedRaff[d] - U3Init[d-1][3])/U3Init[d-1][3])**2 + ((m.UofExt[d] - U1Init[d-1][6])/U1Init[d-1][6])**2 for d in m.Data) + (10**m.rho[1])*(sum(((L2b[i-1] - m.b[i])/L2b[i-1])**2 + ((L2K[i-1] - m.K[i])/L2K[i-1])**2 + ((L2Kap[i-1] - m.Kap[i])/L2Kap[i-1])**2 for i in m.Comp) + ((m.LRe - 0.05*Colum*ColL)/(0.05*Colum*ColL))**2)
    else:
        print("ERROR: Confirm Isotype definition")
instance.obj2 = Objective(rule = obj2_expr)

# -------------------------------------------------------------------
# Relax constraints
# -------------------------------------------------------------------

# To decrease the degree of freedom without setting constraints, some constraints are deactivated and some variables are fixed at 0

instance.Ucon.deactivate()

for d in instance.Data:
    instance.UofIC[d].fix(0)
    instance.UofRaff[d].fix(0)
    instance.UofFeedRaff[d].fix(0)
    instance.UofExt[d].fix(0)

for d in instance.Data:
    for i in instance.Comp:
        instance.averageCE[d,i].fix(0)
        instance.averageCR[d,i].fix(0)

# -------------------------------------------------------------------
# Selecting Objective
# -------------------------------------------------------------------

instance.obj1.activate()
instance.obj2.deactivate()

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------
opt = SolverFactory('ipopt')

opt.options['mu_init'] = 1e-4

opt.options['ma57_pivtol'] = 1e-8

opt.options['max_iter'] = 5000

#opt.options['linear_solver'] = 'ma97'
# opt.options['linear_solver'] = 'ma57'
opt.options['linear_solver'] = 'ma27'

# -------------------------------------------------------------------
# Solve successively
# -------------------------------------------------------------------

instance.LRe.fix(0.05*Colum*ColL)

if Isotype == "Henry":
    pass
elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
    instance.b[1].fix(bInit[0]/10.0)
    instance.b[2].fix(bInit[1]/10.0)
else:
    print("ERROR: Confirm Isotype definition")

instance.preprocess()
results = opt.solve(instance, tee=True)
instance.load(results)

if Isotype == "Henry":
    pass
elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
    instance.b[1].free()
    instance.b[2].free()
else:
    print("ERROR: Confirm Isotype definition")

# -------------------------------------------------------------------
# Display Results
# -------------------------------------------------------------------

print("\n-------------------Initial Condition Solution----------------------\n")

for d in instance.Data:
    print("Data%s"%(d))
    for i in instance.Comp:
        print("recovery in Extract %s : \t %s" %(i, value(100*sum(instance.intCE[d,i,j] for j in instance.Col )/sum(instance.UF[d,j,k]*instance.CF[d,i]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

    for i in instance.Comp:
        print("purity in Extract %s : \t %s" %(i, value(100*sum(instance.intCE[d,i,j] for j in instance.Col )/sum(instance.intCE[d,n,j] for n in instance.Comp for j in instance.Col ))))

    for i in instance.Comp:
        print("recovery in Raffinate %s : \t %s" %(i, value(100*sum(instance.intCR[d,i,j] for j in instance.Col )/sum(instance.UF[d,j,k]*instance.CF[d,i]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

    for i in instance.Comp:
        print("purity in Raffinate %s : \t %s" %(i, value(100*sum(instance.intCR[d,i,j] for j in instance.Col )/sum(instance.intCR[d,n,j] for n in instance.Comp for j in instance.Col ))))

for d in instance.Data:
    print("Step Time(Data%s) %s" %(d,value(instance.StepTime[d])))

for d in instance.Data:
    print("Data%s" %(d))
    print("Comp A conc in Extract = %s" %(value((sum(instance.intCE[d,1,j] for j in instance.Col))/(sum(instance.UE[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))
    print("Comp B conc in Extract = %s" %(value((sum(instance.intCE[d,2,j] for j in instance.Col))/(sum(instance.UE[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))
    print("Comp A conc in Raffinate = %s" %(value((sum(instance.intCR[d,1,j] for j in instance.Col))/(sum(instance.UR[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))
    print("Comp B conc in Raffinate = %s" %(value((sum(instance.intCR[d,2,j] for j in instance.Col))/(sum(instance.UR[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))

# -------------------------------------------------------------------
# Setting fitting parameters free
# -------------------------------------------------------------------

for d in instance.Data:
    for i in instance.t:
        instance.U[d,1,i].free()
        instance.U[d,Zone+1,i].free()
        instance.U[d,2*Zone+1,i].free()
        instance.U[d,3*Zone+1,i].free()

for i in instance.Comp:
    instance.K[i].free()
    instance.Kap[i].free()
    if Isotype == "Henry":
        pass
    elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
        instance.b[i].free()
    else:
        print("ERROR: Confirm Isotype definition")

instance.LRe.free()

for d in instance.Data:
    for i in instance.Comp:
        instance.averageCE[d,i].free()
        instance.averageCR[d,i].free()

# -------------------------------------------------------------------
# Selecting Objective
# -------------------------------------------------------------------

instance.obj1.deactivate()
instance.obj2.activate()

# -------------------------------------------------------------------
# Setting flow rate of ASMB
# -------------------------------------------------------------------

for d in instance.Data:
    instance.UofIC[d].free()
    instance.UofRaff[d].free()
    instance.UofFeedRaff[d].free()
    instance.UofExt[d].free()

for d in instance.Data:
    for i in range(0*Zone+1,4*Zone+1):
        # IC0
        instance.Ucon.add(instance.UofIC[d] == instance.U[d,i,instance.t[1]])
        for j in range(1,sum(HTR[0:1])):
            instance.Ucon.add(instance.UofIC[d] == instance.U[d,i,instance.t[1+j*instance.NCP+1]])
        # IC1
        for j in range(sum(HTR[0:3]),sum(HTR[0:4])):
            instance.Ucon.add(instance.UofIC[d] == instance.U[d,i,instance.t[1+j*instance.NCP+1]])
        # IC2
        for j in range(sum(HTR[0:5]),sum(HTR[0:6])):
            instance.Ucon.add(instance.UofIC[d] == instance.U[d,i,instance.t[1+j*instance.NCP+1]])

# Raff
# Desorbent is only supplied from inlet of first colmun in Zone 1
# Raffinate is only collected from outlet of last colmun in Zone 3
for d in instance.Data:
    for j in range(sum(HTR[0:1]),sum(HTR[0:2])):
        for i in range(0*Zone+1,3*Zone+1):
            instance.Ucon.add(instance.UofRaff[d] == instance.U[d,i,instance.t[1+j*instance.NCP+1]])
        for i in range(3*Zone+1,4*Zone+1):
            instance.Ucon.add(instance.Usmall == instance.U[d,i,instance.t[1+j*instance.NCP+1]])

# FeedRaff
# Feed is only supplied from inlet of first colmun in Zone 3
# Raffinate is only collected from outlet of last colmun in Zone 3
for d in instance.Data:
    for j in range(sum(HTR[0:2]),sum(HTR[0:3])):
        for i in range(0*Zone+1,2*Zone+1):
            instance.Ucon.add(instance.Usmall == instance.U[d,i,instance.t[1+j*instance.NCP+1]])
        for i in range(2*Zone+1,3*Zone+1):
            instance.Ucon.add(instance.UofFeedRaff[d] == instance.U[d,i,instance.t[1+j*instance.NCP+1]])
        for i in range(3*Zone+1,4*Zone+1):
            instance.Ucon.add(instance.Usmall == instance.U[d,i,instance.t[1+j*instance.NCP+1]])

# Ext
# Desorbent is only supplied from inlet of first colmun in Zone 1
# Extract is only collected from outlet of last colmun in Zone 1
for d in instance.Data:
    for j in range(sum(HTR[0:4]),sum(HTR[0:5])):
        for i in range(0*Zone+1,1*Zone+1):
            instance.Ucon.add(instance.UofExt[d] == instance.U[d,1,instance.t[1+j*instance.NCP+1]])
        for i in range(1*Zone+1,4*Zone+1):
            instance.Ucon.add(instance.Usmall == instance.U[d,i,instance.t[1+j*instance.NCP+1]])

# -------------------------------------------------------------------
# Caluculating average concentration in each product
# -------------------------------------------------------------------

# Adding constraints to change averageCE and CR from free variables to bound variables

def AveCE_rule(m, d,i):
    return m.averageCE[d,i]*sum(m.UE[d,1,l]*m.HT[d,l]/m.NCP for l in m.t if l != 0) == m.intCE[d,i,1]
instance.AveCE = Constraint(instance.Data,instance.Comp,rule=AveCE_rule)

def AveCR_rule(m, d,i):
    return m.averageCR[d,i]*sum(m.UR[d,3,l]*m.HT[d,l]/m.NCP for l in m.t if l != 0) == m.intCR[d,i,3]
instance.AveCR = Constraint(instance.Data,instance.Comp,rule=AveCR_rule)

# -------------------------------------------------------------------
# Activating constraints
# -------------------------------------------------------------------

instance.Ucon.activate()

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------
opt.options['mu_init'] = 1e-3

opt.options['ma57_pivtol'] = 1e-8

opt.options['max_iter'] = 5000

#opt.options['linear_solver'] = 'ma97'
#opt.options['linear_solver'] = 'ma57'
opt.options['linear_solver'] = 'ma27'

instance.preprocess()
results = opt.solve(instance, tee=True)
instance.load(results)

# -------------------------------------------------------------------
# Display Results
# -------------------------------------------------------------------

print("\n-------------------Estimation Solution----------------------\n")

for d in instance.Data:
    print("Data%s"%(d))
    for i in instance.Comp:
        print("recovery in Extract %s : \t %s" %(i, value(100*sum(instance.intCE[d,i,j] for j in instance.Col )/sum(instance.UF[d,j,k]*instance.CF[d,i]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

    for i in instance.Comp:
        print("purity in Extract %s : \t %s" %(i, value(100*sum(instance.intCE[d,i,j] for j in instance.Col )/sum(instance.intCE[d,n,j] for n in instance.Comp for j in instance.Col ))))

    for i in instance.Comp:
        print("recovery in Raffinate %s : \t %s" %(i, value(100*sum(instance.intCR[d,i,j] for j in instance.Col )/sum(instance.UF[d,j,k]*instance.CF[d,i]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0))))

    for i in instance.Comp:
        print("purity in Raffinate %s : \t %s" %(i, value(100*sum(instance.intCR[d,i,j] for j in instance.Col )/sum(instance.intCR[d,n,j] for n in instance.Comp for j in instance.Col ))))

for d in instance.Data:
    print("Step Time(Data%s) %s" %(d,value(instance.StepTime[d])))

print("\n")

for d in instance.Data:
    print("\nData%s" %(d))
    print("Throughput: %s [m/hr]" %(value(3600*sum(instance.UF[d,j,k]*instance.HT[d,k]/instance.NCP for k in instance.t for j in instance.Col if k != 0))))
    print("Desorbent: %s [m/hr]" %(value(3600*sum(instance.UD[d,j,k]*instance.HT[d,k]/instance.NCP for k in instance.t for j in instance.Col if k != 0))))
    print("D/F ratio: %s" %(value((sum(instance.UD[d,j,k]*instance.HT[d,k]/instance.NCP for k in instance.t for j in instance.Col if k != 0)/(sum(instance.UF[d,j,k]*instance.HT[d,k]/instance.NCP for k in instance.t for j in instance.Col if k != 0))))))
    print("\nComp A conc in extract = %s" %(value((sum(instance.intCE[d,1,j] for j in instance.Col))/(sum(instance.UE[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))
    print("Comp B conc in extract = %s" %(value((sum(instance.intCE[d,2,j] for j in instance.Col))/(sum(instance.UE[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))
    print("Comp A conc in raffinate = %s" %(value((sum(instance.intCR[d,1,j] for j in instance.Col))/(sum(instance.UR[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))
    print("Comp B conc in raffinate = %s" %(value((sum(instance.intCR[d,2,j] for j in instance.Col))/(sum(instance.UR[d,j,k]*instance.HT[d,k]/instance.NCP for j in instance.Col for k in instance.t if k != 0)))))
    print("\nUofIC = %s " %(value(instance.UofIC[d])))
    print("UofRaff = %s " %(value(instance.UofRaff[d])))
    print("UofFeedRaff = %s " %(value(instance.UofFeedRaff[d])))
    print("UofExt = %s " %(value(instance.UofExt[d])))
    print("Fitting term = %s " %(value((instance.CbarE[d,1] - instance.averageCE[d,1])**2 + (instance.CbarR[d,1] - instance.averageCR[d,1])**2 + (instance.CbarE[d,2] - instance.averageCE[d,2])**2 + (instance.CbarR[d,2] - instance.averageCR[d,2])**2)))

print("FlowRateError = (R,Raff,FeedRaff,Ext)")
for d in instance.Data:
    print("Day%s = (%s,%s,%s,%s)" %(d,value((instance.UofIC[d]-UICexp[d-1])*100/UICexp[d-1]),value((instance.UofRaff[d]-URaffexp[d-1])*100/URaffexp[d-1]),value((instance.UofFeedRaff[d]-UFeedRaffexp[d-1])*100/UFeedRaffexp[d-1]),value((instance.UofExt[d]-UExtexp[d-1])*100/UExtexp[d-1])))

print("\n")

print("Recycle Length = %s " %(value(instance.LRe)))
for i in sorted(instance.Comp):
    print("Kap[%s] = %s" %(i,value(instance.Kap[i])))
    print("K[%s] = %s" %(i,value(instance.K[i])))
    if Isotype == "Henry":
        pass
    elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
        print("b[%s] = %s" %(i,value(instance.b[i])))
    else:
        print("ERROR: Confirm Isotype definition")

print("\nEstimation Completed")

# -------------------------------------------------------------------
# Plot graphs
# -------------------------------------------------------------------

dirname01 = dirname00 + f"Plotfile_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
os.makedirs(dirname01, exist_ok=True)

Col = []
t = []
X = []
Comp = []
CA = [[] for i in range(NData)]
CB = [[] for i in range(NData)]
S = np.zeros(NData)

for i in sorted(instance.Col):
	Col.append(i)

for j in range(NData):
    for i in sorted(instance.FET):
        S[j] = S[j] + value(instance.HT[j+1,instance.t[i*instance.NCP+1]])

for i in sorted(instance.Comp):
	Comp.append(i)

for k in range(NData):
    for i in sorted(instance.Col):
        for j in sorted(instance.x):
            CA[k].append(value(instance.C[k+1, 1, i, 1.0, j])/value(instance.CF[k+1,1]))
            CB[k].append(value(instance.C[k+1, 2, i, 1.0, j])/value(instance.CF[k+1,2]))

for i in sorted(instance.Col):
    for j in sorted(instance.x):
        X.append(j + ColL*(i-1))

filename = [[] for i in range(NData)]
for i in range(NData):
    plt.figure(i+1)
    plt.plot(X, CA[i], 'o-', label='Comp A')
    plt.plot(X, CB[i], 'x-', label='Comp B')
    plt.legend()
    plt.xlim(0, ColL*Colum)
    plt.ylim(0, 2)
    plt.xlabel('x [m]')
    plt.ylabel('C/C(feed) [-]')
    plt.title(f"Data{i+1}_Concentration profile")
    filename[i] = dirname01 + f"Data{i+1}_concentration_profile"
    plt.savefig(filename[i], dpi=500)

dirname02 = dirname00 + f"Animationfile_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
os.makedirs(dirname02, exist_ok=True)

for l in range(NData):
    C1all=[]
    C2all=[]

    for k in sorted(instance.t):
        C1 = []
        C2 = []
        for i in sorted(instance.Col):
            for j in sorted(instance.x):
                C1.append((value(instance.C[l+1,1, i, k ,j]))/value(instance.CF[l+1,1]))
                C2.append((value(instance.C[l+1,2, i, k ,j]))/value(instance.CF[l+1,2]))
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
        plt.ylim(0, 2.0)
        plt.xlim(0, ColL*Colum)
        plt.plot(x,u1[i,:], 'o-', label='Comp A')
        plt.plot(x,u2[i,:], 'x-', label='Comp B')
        plt.xlabel('x [m]')
        plt.ylabel('C/C(feed) [-]')
        plt.title(f"Data{l+1} ASMB Concentration profile")
        plt.legend()


    nframe = Nfet*value(instance.NCP)+1
    nframe = nframe*4
    ani = animation.FuncAnimation(fig,update,fargs = (C1_ani,C2_ani,X,'Wave motion'),interval=100,frames=nframe)
    filename = dirname02 + f"Data{l+1}_ASMB_Concentration_Animation.gif"
    ani.save(filename, writer="pillow", dpi=500)