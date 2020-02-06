# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# ASMB model was constructed by Kensuke Suzuki in 2020
# The method used in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

# ==================================================================
# Pyomo script file
# ==================================================================
from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import SolverFactory
from ParamEst import m, data
from ParamEst_Initdata import *
import datetime
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import time

t1 = time.time()

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

Nfex = Nfex+1
Zone = Colum/4

m.Ucon = ConstraintList() # list of constraints on fluid velocity (used in parameter estimation step)

m.intCR = Var(m.Data, m.Comp, m.Col) # integrated concentration of Raffinate
m.intCE = Var(m.Data, m.Comp, m.Col) # integrated concentration of Extract

m.UofR = Var(m.Data) # fluid velocity of Circulation-process in ASMB (used in parameter estimation step)
m.UofRaff = Var(m.Data) # fluid velocity of Raffinate-collected-process in ASMB
m.UofFeedRaff = Var(m.Data) # fluid velocity of Feed-injected-and-Raffinate-collected-process in ASMB
m.UofExt = Var(m.Data) # fluid velocity of Extract-collected-process in ASMB
m.Usmall = Param(initialize = 0.000001) # please reffer to the document. 

instance = m.create_instance(data)


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
#discretizer2.apply_to(instance,nfe=Nfex-1, wrt=instance.x, scheme='BACKWARD')
discretizer2.apply_to(instance,nfe=Nfex-1, wrt=instance.x, scheme='CENTRAL')

discretizer2.apply_to(instance,nfe=Ncstr, wrt=instance.xRe, scheme='BACKWARD')

# -------------------------------------------------------------------
# Constraints of fluid velocity for each collocation point in finite elements
# -------------------------------------------------------------------
for d in sorted(instance.Data):
    for i in range(1,Colum+1):
        instance.Ucon.add(instance.U[d,(i-1)*Zone+1,instance.t[1]] == instance.U[d,(i-1)*Zone+1,instance.t[2]])
        for j in range(1,Nfet+1):
            for k in range(2,value(instance.NCP)+1):
                instance.Ucon.add(instance.U[d,(i-1)*Zone+1,instance.t[(j-1)*instance.NCP+k]] == instance.U[d,(i-1)*Zone+1,instance.t[(j-1)*instance.NCP+k+1]])

# -------------------------------------------------------------------
# Fixing variables using Build Action
# -------------------------------------------------------------------
def Kap_Fix_rule(m):
  for i in sorted(m.Comp):
      m.Kap[i].fix(KapInit[i-1])
instance.Kap_Fix = BuildAction(rule = Kap_Fix_rule)

def K_Fix_rule(m):
  for i in sorted(m.Comp):
      m.K[i].fix(KInit[i-1])
instance.K_Fix = BuildAction(rule = K_Fix_rule)

def b_Fix_rule(m):
  for i in sorted(m.Comp):
      m.b[i].fix(bInit[i-1])
instance.b_Fix = BuildAction(rule = b_Fix_rule)

def U1Fix_rule(m):
    for d in sorted(m.Data):
        m.U[d,1,0].fix(U1Init[d-1][0])
        for i in range(1,m.NCP+1):
            for j in range(0,Nfet):
                m.U[d,1,m.t[m.NCP*j+i+1]].fix(U1Init[d-1][j])
instance.U1Fix = BuildAction(rule = U1Fix_rule)

def U2Fix_rule(m):
    for d in sorted(m.Data):
        m.U[d,2,0].fix(U2Init[d-1][0])
        for i in range(1,m.NCP+1):
            for j in range(0,Nfet):
                m.U[d,Zone+1,m.t[m.NCP*j+i+1]].fix(U2Init[d-1][j])
instance.U2Fix = BuildAction(rule = U2Fix_rule)

def U3Fix_rule(m):
    for d in sorted(m.Data):
        m.U[d,3,0].fix(U3Init[d-1][0])
        for i in range(1,m.NCP+1):
            for j in range(0,Nfet):
                m.U[d,2*Zone+1,m.t[m.NCP*j+i+1]].fix(U3Init[d-1][j])
instance.U3Fix = BuildAction(rule = U3Fix_rule)

def U4Fix_rule(m):
    for d in sorted(m.Data):
        m.U[d,4,0].fix(U4Init[d-1][0])
        for i in range(1,m.NCP+1):
            for j in range(0,Nfet):
                m.U[d,3*Zone+1,m.t[m.NCP*j+i+1]].fix(U4Init[d-1][j])
instance.U4Fix = BuildAction(rule = U4Fix_rule)

def UF_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i!=(2*Zone+1):
                    m.UF[d,i,j].fix(0)
instance.UF_Fix = BuildAction(rule = UF_Fix_rule)

def UE_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i!=Zone:
                    m.UE[d,i,j].fix(0)
instance.UE_Fix = BuildAction(rule = UE_Fix_rule)

def UD_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i!=1:
                    m.UD[d,i,j].fix(0)
instance.UD_Fix = BuildAction(rule = UD_Fix_rule)

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

# Nfet=12,FA:+1
def HT_Fix_rule(m):
    for d in sorted(m.Data):
        for i in m.t:
            if i < value(m.t[m.NCP+2]):
                m.HT[d,i].fix(1.0*HT1[d-1]/StepTimeInit[d-1])
            if value(m.t[m.NCP+1]) < i and i < value(m.t[2*m.NCP+2]):
                m.HT[d,i].fix(1.0*HT2[d-1]/StepTimeInit[d-1])
            if value(m.t[2*m.NCP+1]) < i and i < value(m.t[5*m.NCP+2]):
                m.HT[d,i].fix(1.0*HT3[d-1]/(StepTimeInit[d-1]*3))
            if value(m.t[5*m.NCP+1]) < i and i < value(m.t[6*m.NCP+2]):
                m.HT[d,i].fix(1.0*HT4[d-1]/StepTimeInit[d-1])
            if value(m.t[6*m.NCP+1]) < i and i < value(m.t[7*m.NCP+2]):
                m.HT[d,i].fix(1.0*HT5[d-1]/StepTimeInit[d-1])
            if value(m.t[7*m.NCP+1]) < i:
                m.HT[d,i].fix(1.0*HT6[d-1]/(StepTimeInit[d-1]*5))
instance.HT_Fix = BuildAction(rule = HT_Fix_rule)

# -------------------------------------------------------------------
# Defining calculation of Integrals
# -------------------------------------------------------------------

def _intCR(m, d,i,j):
    return sum(m.UR[d,j,m.t[(k-1)*m.NCP+1+l]]*m.HT[d,m.t[(k-1)*m.NCP+1+l]]*(m.A[l,3]*m.C[d,i,j,m.t[(k-1)*m.NCP+1+l],value(m.L)]) for l in range(1,m.NCP+1) for k in range(1,Nfet+1)) == m.intCR[d,i,j]
instance.intCRrule = Constraint(instance.Data, instance.Comp, instance.Col,rule=_intCR)

def _intCE(m, d,i,j):
    return sum(m.UE[d,j,m.t[(k-1)*m.NCP+1+l]]*m.HT[d,m.t[(k-1)*m.NCP+1+l]]*(m.A[l,3]*m.C[d,i,j,m.t[(k-1)*m.NCP+1+l],value(m.L)]) for l in range(1,m.NCP+1) for k in range(1,Nfet+1)) == m.intCE[d,i,j]
instance.intCErule = Constraint(instance.Data, instance.Comp, instance.Col,rule=_intCE)

# -------------------------------------------------------------------
# Activate if CENTRAL scheme is selected
# -------------------------------------------------------------------

def C_AxialDerivativeConstraintEnd_rule(m,d,i,j,k):
	return m.dCdx[d,i,j,k,value(m.L)]==(3*m.C[d,i,j,k,value(m.L)]-4*m.C[d,i,j,k,m.x[Nfex-1]]+m.C[d,i,j,k,m.x[Nfex-2]])/(2*m.L/(Nfex-1)) 
instance.C_AxialDerivativeConstraintEnd=Constraint(instance.Data, instance.Comp, instance.Col, instance.t,  rule=C_AxialDerivativeConstraintEnd_rule)

# -------------------------------------------------------------------
# Objective Function
# -------------------------------------------------------------------

def obj1_expr(m):
    return 3600*sum(sum(m.UF[1,j,k]*m.HT[1,k]/m.NCP for k in m.t if k != 0) for j in m.Col)
instance.obj1 = Objective(rule = obj1_expr, sense = maximize)

def obj2_expr(m):
    return (1/len(m.Data))*sum((m.CbarE[d,i] - m.averageCE[d,i])**2 + (m.CbarR[d,i] - m.averageCR[d,i])**2 for d in m.Data for i in m.Comp) + ((10**m.rho[2])/len(m.Data))*sum(((m.UofR[d] - U1Init[d-1][0])/U1Init[d-1][0])**2 + ((m.UofRaff[d] - U1Init[d-1][1])/U1Init[d-1][1])**2 + ((m.UofFeedRaff[d] - U3Init[d-1][3])/U3Init[d-1][3])**2 + ((m.UofExt[d] - U1Init[d-1][6])/U1Init[d-1][6])**2 for d in m.Data) + (10**m.rho[1])*(sum(((L2b[i-1] - m.b[i])/L2b[i-1])**2 + ((L2K[i-1] - m.K[i])/L2K[i-1])**2 + ((L2Kap[i-1] - m.Kap[i])/L2Kap[i-1])**2 for i in m.Comp) + ((m.LRe - 0.05*Colum*ColL)/(0.05*Colum*ColL))**2) 
instance.obj2 = Objective(rule = obj2_expr)

# -------------------------------------------------------------------
# Relax constraints
# -------------------------------------------------------------------

instance.Ucon.deactivate()

for d in instance.Data:
    instance.UofR[d].fix(0)
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
# Solver Options
# -------------------------------------------------------------------

# Create the ipopt_sens solver plugin using the ASL interface
solver = 'ipopt_sens'
solver_io = 'nl'
stream_solver = True    # True prints solver output to screen
keepfiles =     False    # True prints intermediate file names (.nl,.sol,...)
opt = SolverFactory(solver,solver_io=solver_io)

opt.options['mu_init'] = 1e-4
           
opt.options['ma57_pivtol'] = 1e-8

opt.options['max_iter'] = 5000

# -------------------------------------------------------------------
# Solve successively
# -------------------------------------------------------------------

instance.LRe.fix(0.05*Colum*ColL)

instance.b[1].fix(bInit[0]/10.0)
instance.b[2].fix(bInit[1]/10.0)

instance.preprocess()
results = opt.solve(instance, tee=True)
instance.load(results)

instance.b[1].free()
instance.b[2].free()

# -------------------------------------------------------------------
# Display Results
# -------------------------------------------------------------------

# print("\n-------------------Initial Condition Solution----------------------\n")

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
    instance.b[i].free()

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
    instance.UofR[d].free()
    instance.UofRaff[d].free()
    instance.UofFeedRaff[d].free()
    instance.UofExt[d].free()

# R0,R1,R2
for d in instance.Data:
    for i in range(1,Colum+1):
        instance.Ucon.add(instance.UofR[d] == instance.U[d,i,instance.t[1]])
        instance.Ucon.add(instance.UofR[d] == instance.U[d,i,instance.t[5*instance.NCP+2]])
        for j in range(7,Nfet):
            instance.Ucon.add(instance.UofR[d] == instance.U[d,i,instance.t[j*instance.NCP+2]])

# Raff
for d in instance.Data:
    for i in range(1,Colum):
        instance.Ucon.add(instance.UofRaff[d] == instance.U[d,i,instance.t[instance.NCP+2]])
    instance.Ucon.add(instance.Usmall == instance.U[d,4,instance.t[instance.NCP+2]])

# FeedRaff
for d in instance.Data:
    for j in range(2,5):
        instance.Ucon.add(instance.UofFeedRaff[d] == instance.U[d,2*Zone+1,instance.t[j*instance.NCP+2]])
        instance.Ucon.add(instance.Usmall == instance.U[d,1,instance.t[j*instance.NCP+2]])
        instance.Ucon.add(instance.Usmall == instance.U[d,2,instance.t[j*instance.NCP+2]])
        instance.Ucon.add(instance.Usmall == instance.U[d,4,instance.t[j*instance.NCP+2]])

# Ext
for d in instance.Data:
    instance.Ucon.add(instance.UofExt[d] == instance.U[d,Zone,instance.t[6*instance.NCP+2]])
    for i in range(2,Colum+1):
        instance.Ucon.add(instance.Usmall == instance.U[d,i,instance.t[6*instance.NCP+2]])

# -------------------------------------------------------------------
# Caluculating average concentration
# -------------------------------------------------------------------

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

opt.options['ma57_pivtol'] = 1e-7

opt.options['max_iter'] = 5000

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------

instance.eta1 = Var()

nominal_eta1   = 5
perturbed_eta1 = 5

instance.consteta1 = Constraint(expr=instance.eta1 == nominal_eta1)

# declare suffixes
instance.sens_state_0 = Suffix(direction=Suffix.EXPORT)
instance.sens_state_1 = Suffix(direction=Suffix.EXPORT)
instance.sens_state_value_1 = Suffix(direction=Suffix.EXPORT)
instance.sens_sol_state_1  = Suffix(direction=Suffix.IMPORT)
instance.sens_init_constr  = Suffix(direction=Suffix.EXPORT)
instance.red_hessian = Suffix(direction=Suffix.EXPORT)

# set sIPOPT data
opt.options['run_sens'] = 'yes'
opt.options['compute_red_hessian'] = 'yes'

instance.sens_state_0[instance.eta1] = 1
instance.sens_state_1[instance.eta1] = 1
instance.sens_state_value_1[instance.eta1] = perturbed_eta1
instance.sens_init_constr[instance.consteta1] = 1

instance.red_hessian[instance.Kap[1]] = 1
instance.red_hessian[instance.K[1]] = 2
instance.red_hessian[instance.b[1]] = 3
instance.red_hessian[instance.Kap[2]] = 4
instance.red_hessian[instance.K[2]] = 5
instance.red_hessian[instance.b[2]] = 6
instance.red_hessian[instance.LRe] = 7
instance.red_hessian[instance.UofR[1]] = 8
instance.red_hessian[instance.UofRaff[1]] = 9
instance.red_hessian[instance.UofFeedRaff[1]] = 10
instance.red_hessian[instance.UofExt[1]] = 11
instance.red_hessian[instance.UofR[2]] = 12
instance.red_hessian[instance.UofRaff[2]] = 13
instance.red_hessian[instance.UofFeedRaff[2]] = 14
instance.red_hessian[instance.UofExt[2]] = 15
instance.red_hessian[instance.UofR[3]] = 16
instance.red_hessian[instance.UofRaff[3]] = 17
instance.red_hessian[instance.UofFeedRaff[3]] = 18
instance.red_hessian[instance.UofExt[3]] = 19
instance.red_hessian[instance.UofR[4]] = 20
instance.red_hessian[instance.UofRaff[4]] = 21
instance.red_hessian[instance.UofFeedRaff[4]] = 22
instance.red_hessian[instance.UofExt[4]] = 23
instance.red_hessian[instance.UofR[5]] = 24
instance.red_hessian[instance.UofRaff[5]] = 25
instance.red_hessian[instance.UofFeedRaff[5]] = 26
instance.red_hessian[instance.UofExt[5]] = 27
instance.red_hessian[instance.UofR[6]] = 28
instance.red_hessian[instance.UofRaff[6]] = 29
instance.red_hessian[instance.UofFeedRaff[6]] = 30
instance.red_hessian[instance.UofExt[6]] = 31
instance.red_hessian[instance.UofR[7]] = 32
instance.red_hessian[instance.UofRaff[7]] = 33
instance.red_hessian[instance.UofFeedRaff[7]] = 34
instance.red_hessian[instance.UofExt[7]] = 35
instance.red_hessian[instance.UofR[8]] = 36
instance.red_hessian[instance.UofRaff[8]] = 37
instance.red_hessian[instance.UofFeedRaff[8]] = 38
instance.red_hessian[instance.UofExt[8]] = 39
instance.red_hessian[instance.UofR[9]] = 40
instance.red_hessian[instance.UofRaff[9]] = 41
instance.red_hessian[instance.UofFeedRaff[9]] = 42
instance.red_hessian[instance.UofExt[9]] = 43
instance.red_hessian[instance.UofR[10]] = 44
instance.red_hessian[instance.UofRaff[10]] = 45
instance.red_hessian[instance.UofFeedRaff[10]] = 46
instance.red_hessian[instance.UofExt[10]] = 47
instance.red_hessian[instance.UofR[11]] = 48
instance.red_hessian[instance.UofRaff[11]] = 49
instance.red_hessian[instance.UofFeedRaff[11]] = 50
instance.red_hessian[instance.UofExt[11]] = 51
instance.red_hessian[instance.UofR[12]] = 52
instance.red_hessian[instance.UofRaff[12]] = 53
instance.red_hessian[instance.UofFeedRaff[12]] = 54
instance.red_hessian[instance.UofExt[12]] = 55
instance.red_hessian[instance.UofR[13]] = 56
instance.red_hessian[instance.UofRaff[13]] = 57
instance.red_hessian[instance.UofFeedRaff[13]] = 58
instance.red_hessian[instance.UofExt[13]] = 59
instance.red_hessian[instance.UofR[14]] = 60
instance.red_hessian[instance.UofRaff[14]] = 61
instance.red_hessian[instance.UofFeedRaff[14]] = 62
instance.red_hessian[instance.UofExt[14]] = 63
instance.red_hessian[instance.UofR[15]] = 64
instance.red_hessian[instance.UofRaff[15]] = 65
instance.red_hessian[instance.UofFeedRaff[15]] = 66
instance.red_hessian[instance.UofExt[15]] = 67
instance.red_hessian[instance.UofR[16]] = 68
instance.red_hessian[instance.UofRaff[16]] = 69
instance.red_hessian[instance.UofFeedRaff[16]] = 70
instance.red_hessian[instance.UofExt[16]] = 71
instance.red_hessian[instance.UofR[17]] = 72
instance.red_hessian[instance.UofRaff[17]] = 73
instance.red_hessian[instance.UofFeedRaff[17]] = 74
instance.red_hessian[instance.UofExt[17]] = 75
instance.red_hessian[instance.UofR[18]] = 76
instance.red_hessian[instance.UofRaff[18]] = 77
instance.red_hessian[instance.UofFeedRaff[18]] = 78
instance.red_hessian[instance.UofExt[18]] = 79
instance.red_hessian[instance.UofR[19]] = 80
instance.red_hessian[instance.UofRaff[19]] = 81
instance.red_hessian[instance.UofFeedRaff[19]] = 82
instance.red_hessian[instance.UofExt[19]] = 83
instance.red_hessian[instance.UofR[20]] =	84
instance.red_hessian[instance.UofRaff[20]] = 85
instance.red_hessian[instance.UofFeedRaff[20]] = 86
instance.red_hessian[instance.UofExt[20]] = 87
instance.red_hessian[instance.UofR[21]] =	 88
instance.red_hessian[instance.UofRaff[21]] = 89
instance.red_hessian[instance.UofFeedRaff[21]] = 90
instance.red_hessian[instance.UofExt[21]] = 91
instance.red_hessian[instance.UofR[22]] =	 92
instance.red_hessian[instance.UofRaff[22]] = 93
instance.red_hessian[instance.UofFeedRaff[22]] = 94
instance.red_hessian[instance.UofExt[22]] = 95
instance.red_hessian[instance.UofR[23]] =	 96
instance.red_hessian[instance.UofRaff[23]] = 97
instance.red_hessian[instance.UofFeedRaff[23]] = 98
instance.red_hessian[instance.UofExt[23]] = 99
instance.red_hessian[instance.UofR[24]] =	 100
instance.red_hessian[instance.UofRaff[24]] = 101
instance.red_hessian[instance.UofFeedRaff[24]] = 102
instance.red_hessian[instance.UofExt[24]] = 103
instance.red_hessian[instance.UofR[25]] =	 104
instance.red_hessian[instance.UofRaff[25]] = 105
instance.red_hessian[instance.UofFeedRaff[25]] = 106
instance.red_hessian[instance.UofExt[25]] = 107
instance.red_hessian[instance.UofR[26]] =	 108
instance.red_hessian[instance.UofRaff[26]] = 109
instance.red_hessian[instance.UofFeedRaff[26]] = 110
instance.red_hessian[instance.UofExt[26]] = 111
instance.red_hessian[instance.UofR[27]] =	 112
instance.red_hessian[instance.UofRaff[27]] = 113
instance.red_hessian[instance.UofFeedRaff[27]] = 114
instance.red_hessian[instance.UofExt[27]] = 115
instance.red_hessian[instance.UofR[28]] =	 116
instance.red_hessian[instance.UofRaff[28]] = 117
instance.red_hessian[instance.UofFeedRaff[28]] = 118
instance.red_hessian[instance.UofExt[28]] = 119
instance.red_hessian[instance.UofR[29]] =	 120
instance.red_hessian[instance.UofRaff[29]] = 121
instance.red_hessian[instance.UofFeedRaff[29]] = 122
instance.red_hessian[instance.UofExt[29]] = 123
instance.red_hessian[instance.UofR[30]] =	 124
instance.red_hessian[instance.UofRaff[30]] = 125
instance.red_hessian[instance.UofFeedRaff[30]] = 126
instance.red_hessian[instance.UofExt[30]] = 127

# Send the model to ipopt_sens and collect the solution
instance.preprocess()
results = opt.solve(instance, keepfiles=keepfiles, tee=stream_solver)
instance.load(results)

# =============================================================================
# 
# =============================================================================

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
    print("\nUofR = %s " %(value(instance.UofR[d])))
    print("UofRaff = %s " %(value(instance.UofRaff[d])))
    print("UofFeedRaff = %s " %(value(instance.UofFeedRaff[d])))
    print("UofExt = %s " %(value(instance.UofExt[d])))
    print("Fitting term = %s " %(value((instance.CbarE[d,1] - instance.averageCE[d,1])**2 + (instance.CbarR[d,1] - instance.averageCR[d,1])**2 + (instance.CbarE[d,2] - instance.averageCE[d,2])**2 + (instance.CbarR[d,2] - instance.averageCR[d,2])**2)))

print("FlowRateError = (R,Raff,FeedRaff,Ext)")
for d in instance.Data:
    print("Day%s = (%s,%s,%s,%s)" %(d,value((instance.UofR[d]-URexp[d-1])*100/URexp[d-1]),value((instance.UofRaff[d]-URaffexp[d-1])*100/URaffexp[d-1]),value((instance.UofFeedRaff[d]-UFeedRaffexp[d-1])*100/UFeedRaffexp[d-1]),value((instance.UofExt[d]-UExtexp[d-1])*100/UExtexp[d-1])))

print("\n")

print("Recycle Length = %s " %(value(instance.LRe)))
for i in sorted(instance.Comp):
    print("Kap[%s] = %s" %(i,value(instance.Kap[i])))
    print("K[%s] = %s" %(i,value(instance.K[i])))
    print("b[%s] = %s" %(i,value(instance.b[i])))

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
CFru = [[] for i in range(NData)]
CGlu = [[] for i in range(NData)]
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
            CFru[k].append(value(instance.C[k+1, 1, i, 1.0, j])/value(instance.CF[k+1,1]))
            CGlu[k].append(value(instance.C[k+1, 2, i, 1.0, j])/value(instance.CF[k+1,2]))

for i in sorted(instance.Col):
    for j in sorted(instance.x):
        X.append(j + ColL*(i-1))

filename = [[] for i in range(NData)]
for i in range(NData):
    plt.figure(i+1)
    line1, = plt.plot(X, CFru[i], 'bx-', label='Comp A')
    line2, = plt.plot(X, CGlu[i], 'r+-', label='Comp B')
    plt.legend(handles=[line1, line2])
    plt.xlim(0, ColL*Colum)
    plt.ylim(0, 2)
    plt.xlabel('x (m)')
    plt.ylabel('C/C(feed)')
    plt.title(f"Data{i+1}_Concentration profile")
    filename[i] = dirname01 + f"Data{i+1}_concentration_profile"
    plt.savefig(filename[i])
