# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kensuke Suzuki in 2020
# The method used in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.
#                          : Zavala V, Biegler LT, Large-scale parameter estimation in low-density polyethylene tubular reactors. Industrial and Engineering Chemistry Research. 2006;45,23;7867-7881
#                          : Tie et al., Experimental evaluation of simulated moving bed reactor for transesterification reaction synthesis of glycol ether ester. Adsorption. 2019;25,4;795-807

# ==================================================================
# Pyomo script file
# ==================================================================

from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import SolverFactory
from SMB_Model_ParamEst import m
from Initdata_SMB_ParamEst import *
from GaussRadauQuadrature import lglnodes
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import datetime
import sys
import os

dirname00 = f"Prameter_Estimation_Result_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
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

# -------------------------------------------------------------------
# Display implementation details
# -------------------------------------------------------------------

print(f"\n----- Discritization Scheme for Axial Coordinate: {DScheme} -----")
print(f"\n----- Equilibrium Isotherm: {Isotype} -----")
print(f"\n----- LDF model based on {Phase} phase -----")

if Axial_D == True:
    print(f"\n----- Axial dispersion is implemented -----")
elif Axial_D == False:
	pass
else:
    sys.exit("ERROR: Confirm Axial_D definition")

if DV == True:
    print("\n----- Dead Volume is implemented -----")
elif DV == False:
    pass
else:
    sys.exit("ERROR: Confirm DeadVolume definition")

if Tikhonov == True:
    if EVM == True:
        print("\n----- Tikhonov regularization and EVM formulation are implemented -----")
    elif EVM == False:
        print("\n----- Tikhonov regularization is implemented -----")
    else:
        sys.exit("ERROR: Confirm EVM definition")
    
elif Tikhonov == False:
    print("\n----- Least square method is implemented -----")
else:
    sys.exit("ERROR: Confirm Tikhonov definition")

# -------------------------------------------------------------------
# Zone Difinition
# -------------------------------------------------------------------

Nfex = Nfex+1
Zone = int(Colum/NZone)

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

discretizer2.apply_to(instance,nfe=Nfex-1, wrt=instance.x, scheme=DScheme)

if DV == True:
    discretizer2.apply_to(instance,nfe=Ncstr, wrt=instance.xRe, scheme='BACKWARD')

# -------------------------------------------------------------------
# Objective Function
# -------------------------------------------------------------------

def obj1_expr(m):
    return sum(m.e[d] for d in m.Data)
instance.obj1 = Objective(rule = obj1_expr, sense = minimize)

# Objective function 2 to estimate model parameters together with fluid velocities

def obj2_expr(m):
    LS_term = (1/len(m.Data))*sum((m.CbarE[d,i] - m.averageCE[d,i])**2 + (m.CbarR[d,i] - m.averageCR[d,i])**2 for d in m.Data for i in m.Comp)
    Tikhonov_term = 0
    EVM_term = 0
    if Tikhonov == True:
        Tikhonov_term += sum(((L2H[i-1] - m.H[i])/L2H[i-1])**2 + ((L2Kap[i-1] - m.Kap[i])/L2Kap[i-1])**2 for i in m.Comp)
        if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
            Tikhonov_term += sum(((L2b[i-1] - m.b[i])/L2b[i-1])**2 for i in m.Comp)
        if Axial_D == True:
            Tikhonov_term += sum(((L2Dax[i-1] - m.Dax[i])/L2Dax[i-1])**2 for i in m.Comp)
        if DV == True:
            Tikhonov_term += ((L2LRe - m.LRe)/L2LRe)**2
        if EVM == True:
            EVM_term += (1/len(m.Data))*sum(((m.UD[d,0*Zone+1,m.t[2]] - UDinit[d-1])/UDinit[d-1])**2 + ((m.UE[d,1*Zone,m.t[2]] - UEinit[d-1])/UEinit[d-1])**2 + 
                                            ((m.UF[d,2*Zone+1,m.t[2]] - UFinit[d-1])/UFinit[d-1])**2 + ((m.U[d,3*Zone+1,m.t[2]] - UReinit[d-1])/UReinit[d-1])**2 for d in m.Data)
    return LS_term + m.rho[1]*Tikhonov_term + m.rho[2]*EVM_term
    
instance.obj2 = Objective(rule = obj2_expr)

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
    if Axial_D == True:
        def C_Axial2ndDerivativeConstraintEnd_rule(m,d,i,j,k):
            return m.dCdx2[d,i,j,k,value(m.L)]==(2*m.C[d,i,j,k,value(m.L)]-5*m.C[d,i,j,k,m.x[Nfex-1]]+4*m.C[d,i,j,k,m.x[Nfex-2]]-m.C[d,i,j,k,m.x[Nfex-3]])/((m.L/(Nfex-1))**3) 
        instance.C_Axial2ndDerivativeConstraintEnd=Constraint(instance.Data, instance.Comp, instance.Col, instance.t,  rule=C_Axial2ndDerivativeConstraintEnd_rule)

# -------------------------------------------------------------------
# Fixing variables using Build Action
# -------------------------------------------------------------------

# ------------------------------------------------------
# Fixing model parameter at initial value
# ------------------------------------------------------

def Kap_Fix_rule(m):
  for i in sorted(m.Comp):
      m.Kap[i].fix(Kapinit[i-1])
instance.Kap_Fix = BuildAction(rule = Kap_Fix_rule)

def H_Fix_rule(m):
  for i in sorted(m.Comp):
      m.H[i].fix(Hinit[i-1])
instance.H_Fix = BuildAction(rule = H_Fix_rule)

if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
    def b_Fix_rule(m):
        for i in sorted(m.Comp):
            m.b[i].fix(binit[i-1])
    instance.b_Fix = BuildAction(rule = b_Fix_rule)

if Axial_D == True:
    def Dax_Fix_rule(m):
        for i in sorted(m.Comp):
            m.Dax[i].fix(Daxinit[i-1])
    instance.Dax_Fix = BuildAction(rule = Dax_Fix_rule)

# ------------------------------------------------------
# Fixing manipulated velocities 
# ------------------------------------------------------

# Feed is only supplied from inlet of first colmun in Zone 3
def UF_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i==(2*Zone+1):
                    m.UF[d,i,j].fix(UFinit[d-1])
                else:
                    m.UF[d,i,j].fix(0)
instance.UF_Fix = BuildAction(rule = UF_Fix_rule)

# Extract is only collected from outlet of last colmun in Zone 1
def UE_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i==(1*Zone):
                    m.UE[d,i,j].fix(UEinit[d-1])
                else:
                    m.UE[d,i,j].fix(0)
instance.UE_Fix = BuildAction(rule = UE_Fix_rule)

# Desorbent is only supplied from inlet of first colmun in Zone 1
def UD_Fix_rule(m):
    for d in sorted(m.Data):
        for j in m.t:
            for i in m.Col:
                if i==(0*Zone+1):
                    m.UD[d,i,j].fix(UDinit[d-1])
                else:
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

def U4Fix_rule(m):
    for d in sorted(m.Data):
        for x in range(3*Zone+1, 4*Zone+1):
            for j in m.t:
                m.U[d,x,j].fix(UReinit[d-1])
instance.U4Fix = BuildAction(rule = U4Fix_rule)

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
            m.HT[d,i].fix(1/Nfet)
instance.HT_Fix = BuildAction(rule = HT_Fix_rule)

# -------------------------------------------------------------------
# Activate and Deactivate Constraints
# -------------------------------------------------------------------

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

#### Create the ipopt_sens solver plugin using the ASL interface
solver = 'ipopt_sens'
solver_io = 'nl'
stream_solver = True    # True prints solver output to screen
keepfiles =     False    # True prints intermediate file names (.nl,.sol,...)
opt = SolverFactory(solver,solver_io=solver_io)
####

opt.options['mu_init'] = 1e-5
#opt.options['ma57_pivtol'] = 1e-8
opt.options['max_iter'] = 5000
# opt.options['linear_system_scaling'] = 'mc19'
#
# opt.options['linear_solver'] = 'ma97'
# opt.options['linear_solver'] = 'ma57'
opt.options['linear_solver'] = 'ma27'

# -------------------------------------------------------------------
# Solve successively
# -------------------------------------------------------------------

if DV == True:
    instance.LRe.fix(LReinit)

for i in range(ite_e):
    for d in instance.Data:
        instance.e[d].fix(1.0*10**-(8-ite_e+i+1))

    print(f"e = 1.0e{-(8-ite_e+i+1)}")

    instance.preprocess()
    results = opt.solve(instance, tee=True)
    instance.load(results)

    instance.e.free()

instance.preprocess()
results = opt.solve(instance, keepfiles=keepfiles, tee=stream_solver)
instance.load(results)

# -------------------------------------------------------------------
# Display Results
# -------------------------------------------------------------------

print("\n-------------------Initial Condition Solution----------------------\n")

print("--------------------------------------------------")
print("Initial parameters")
if DV == True:
    print(f"Recycle Length = {LReinit}")
for i in sorted(instance.Comp):
    print("Kap[%s] = %s" %(i,Kapinit[i-1]))
    print("H[%s] = %s" %(i,Hinit[i-1]))
    if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
        print("b[%s] = %s" %(i,binit[i-1]))
    if Axial_D == True:
        print(f"Dax[{i}] = {Daxinit[i-1]}")
print("--------------------------------------------------")

print("\n")

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

print("\n")

for d in instance.Data:
    print(f"Data {d}: Objective function value (CSS constraint) = {value(instance.e[d])}")

print("\n")

if sum((value(instance.e[d])) for d in instance.Data) >= 3.0e-8:
    sys.exit("ERROR: Simulation not converged")
else:
    print("Simulation converged")

print("\n")

# -------------------------------------------------------------------
# Plot initialized graphs
# -------------------------------------------------------------------
if ini_plots == True:
    dirname0ini = dirname00 + f"InitialPlotfile_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}/"
    os.makedirs(dirname0ini, exist_ok=True)

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
        plt.ylim(0, 1.2)
        plt.xlabel('x [m]')
        plt.ylabel('C/C(feed) [-]')
        plt.title(f"Data{i+1}_Concentration profile")
        filename[i] = dirname0ini + f"Data{i+1}_concentration_profile"
        plt.savefig(filename[i], dpi=500)
elif ini_plots == False:
    pass
else:
    sys.exit("ERROR: Confirm ini_plots definition")

# -------------------------------------------------------------------
# Setting fitting parameters free
# -------------------------------------------------------------------

for d in instance.Data:
    for i in instance.t:
        instance.UF[d,2*Zone+1,i].free()
        instance.UE[d,1*Zone,i].free() 
        instance.UD[d,0*Zone+1,i].free()
        for x in range(3*Zone+1, 4*Zone+1):
            instance.U[d,x,i].free()

for i in instance.Comp:
    instance.H[i].free()
    instance.Kap[i].free()
    if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
        instance.b[i].free()     
    if Axial_D == True:
        instance.Dax[i].free()

for d in instance.Data:
    instance.e[d].fix(0)

if DV == True:
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
# Setting flow rate
# -------------------------------------------------------------------

for d in instance.Data:
    for i in range(1,Nfet):
        instance.U_constraints.add(instance.UF[d,2*Zone+1,instance.t[2]] == instance.UF[d,2*Zone+1,instance.t[i*instance.NCP+2]])
        instance.U_constraints.add(instance.UE[d,1*Zone,instance.t[2]] == instance.UE[d,1*Zone,instance.t[i*instance.NCP+2]])
        instance.U_constraints.add(instance.UD[d,0*Zone+1,instance.t[2]] == instance.UD[d,0*Zone+1,instance.t[i*instance.NCP+2]])
        for x in range(3*Zone+1, 4*Zone+1):
            instance.U_constraints.add(instance.U[d,x,instance.t[2]] == instance.U[d,x,instance.t[i*instance.NCP+2]])

# if HT_Const == True:
#     for d in instance.Data:
#         for i in range(1,Nfet):
#             instance.HT_constraints.add(instance.HT[d,instance.t[2]] == instance.HT[d,instance.t[i*instance.NCP+2]])
# elif HT_Const == False:
#     pass
# else:
#     print("ERROR: Confirm HT_Const definition")

for d in instance.Data:
    instance.StepTime[d].fix(StepTimeInit[d-1])
    for i in range(0,Nfet):
        for j in range(2,value(instance.NCP)+1):
            instance.U_constraints.add(instance.UF[d,2*Zone+1,instance.t[i*instance.NCP+j]] == instance.UF[d,2*Zone+1,instance.t[i*instance.NCP+j+1]])
            instance.U_constraints.add(instance.UE[d,1*Zone,instance.t[i*instance.NCP+j]] == instance.UE[d,1*Zone,instance.t[i*instance.NCP+j+1]])
            instance.U_constraints.add(instance.UD[d,0*Zone+1,instance.t[i*instance.NCP+j]] == instance.UD[d,0*Zone+1,instance.t[i*instance.NCP+j+1]])
            for x in range(3*Zone+1, 4*Zone+1):
                instance.U_constraints.add(instance.U[d,x,instance.t[i*instance.NCP+j]] == instance.U[d,x,instance.t[i*instance.NCP+j+1]])
    instance.U_constraints.add(instance.UF[d,2*Zone+1,instance.t[1]] == instance.UF[d,2*Zone+1,instance.t[2]])
    instance.U_constraints.add(instance.UE[d,1*Zone,instance.t[1]] == instance.UE[d,1*Zone,instance.t[2]])
    instance.U_constraints.add(instance.UD[d,0*Zone+1,instance.t[1]] == instance.UD[d,0*Zone+1,instance.t[2]])
    for x in range(3*Zone+1, 4*Zone+1):
        instance.U_constraints.add(instance.U[d,x,instance.t[1]] == instance.U[d,x,instance.t[2]])

if Tikhonov == True:
    if EVM == False:
        for d in instance.Data:
            instance.UF[d,2*Zone+1,instance.t[2]].fix(UFinit[d-1])
            instance.UE[d,1*Zone,instance.t[2]].fix(UEinit[d-1])
            instance.UD[d,0*Zone+1,instance.t[2]].fix(UDinit[d-1])
            for x in range(3*Zone+1, 4*Zone+1):
                instance.U[d,x,instance.t[2]].fix(UReinit[d-1])
elif Tikhonov == False:
    for d in instance.Data:
        instance.UF[d,2*Zone+1,instance.t[2]].fix(UFinit[d-1])
        instance.UE[d,1*Zone,instance.t[2]].fix(UEinit[d-1])
        instance.UD[d,0*Zone+1,instance.t[2]].fix(UDinit[d-1])
        for x in range(3*Zone+1, 4*Zone+1):
            instance.U[d,x,instance.t[2]].fix(UReinit[d-1])

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

# instance.HTSum.activate()
instance.U_constraints.activate()
# instance.HT_constraints.activate()

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------

#### Create the ipopt_sens solver plugin using the ASL interface
solver = 'ipopt_sens'
solver_io = 'nl'
stream_solver = True    # True prints solver output to screen
keepfiles =     False    # True prints intermediate file names (.nl,.sol,...)
opt = SolverFactory(solver,solver_io=solver_io)
####

opt.options['mu_init'] = 1e-3
#opt.options['ma57_pivtol'] = 1e-8
opt.options['max_iter'] = 5000
# opt.options['linear_system_scaling'] = 'mc19'
#
# opt.options['linear_solver'] = 'ma97'
# opt.options['linear_solver'] = 'ma57'
opt.options['linear_solver'] = 'ma27'

# -------------------------------------------------------------------
# Solver options
# -------------------------------------------------------------------

###
instance.eta1 = Var()

nominal_eta1   = 5
perturbed_eta1 = 5

instance.consteta1 = Constraint(expr=instance.eta1 == nominal_eta1)
###

### declare suffixes
instance.sens_state_0 = Suffix(direction=Suffix.EXPORT)
instance.sens_state_1 = Suffix(direction=Suffix.EXPORT)
instance.sens_state_value_1 = Suffix(direction=Suffix.EXPORT)
instance.sens_sol_state_1  = Suffix(direction=Suffix.IMPORT)
instance.sens_init_constr  = Suffix(direction=Suffix.EXPORT)
instance.red_hessian = Suffix(direction=Suffix.EXPORT)
###

### set sIPOPT data
opt.options['run_sens'] = 'yes'
opt.options['compute_red_hessian'] = 'yes'
###

###
instance.sens_state_0[instance.eta1] = 1
instance.sens_state_1[instance.eta1] = 1
instance.sens_state_value_1[instance.eta1] = perturbed_eta1
instance.sens_init_constr[instance.consteta1] = 1
###

instance.red_hessian[instance.H[1]] = 1
instance.red_hessian[instance.H[2]] = 2
instance.red_hessian[instance.Kap[1]] = 3
instance.red_hessian[instance.Kap[2]] = 4

if Isotype == "Henry":
    if Axial_D == True:
        instance.red_hessian[instance.Dax[1]] = 5
        instance.red_hessian[instance.Dax[2]] = 6
        if DV == True:
            instance.red_hessian[instance.LRe] = 7
    elif Axial_D == False:
        if DV == True:
            instance.red_hessian[instance.LRe] = 5

elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
    instance.red_hessian[instance.b[1]] = 5
    instance.red_hessian[instance.b[2]] = 6
    if Axial_D == True:
        instance.red_hessian[instance.Dax[1]] = 7
        instance.red_hessian[instance.Dax[2]] = 8
        if DV == True:
            instance.red_hessian[instance.LRe] = 9
    elif Axial_D == False:
        if DV == True:
            instance.red_hessian[instance.LRe] = 7


### Send the model to ipopt_sens and collect the solution
instance.preprocess()
results = opt.solve(instance, keepfiles=keepfiles, tee=stream_solver)
instance.load(results)
###

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
    print("\nUF = %s " %(value(instance.UF[d,2*Zone+1,instance.t[2]])))
    print("UE = %s " %(value(instance.UE[d,1*Zone,instance.t[2]])))
    print("UD = %s " %(value(instance.UD[d,0*Zone+1,instance.t[2]])))
    print("URe = %s " %(value(instance.U[d,3*Zone+1,instance.t[2]])))
    print("Fitting term = %s " %(value((instance.CbarE[d,1] - instance.averageCE[d,1])**2 + (instance.CbarR[d,1] - instance.averageCR[d,1])**2 + (instance.CbarE[d,2] - instance.averageCE[d,2])**2 + (instance.CbarR[d,2] - instance.averageCR[d,2])**2)))

print("\n")

if EVM == True:
    print("FlowRateError [%] = (UF,UE,UD,URe)")
    for d in instance.Data:
        print("Day%s = (%s,%s,%s,%s)" %(d,value((instance.UF[d,2*Zone+1,instance.t[2]]-UFinit[d-1])*100/UFinit[d-1]),value((instance.UE[d,1*Zone,instance.t[2]]-UEinit[d-1])*100/UEinit[d-1]),value((instance.UD[d,0*Zone+1,instance.t[2]]-UDinit[d-1])*100/UDinit[d-1]),value((instance.U[d,3*Zone+1,instance.t[2]]-UReinit[d-1])*100/UReinit[d-1])))

print("\n")

print("--------------------------------------------------")
if DV == True:
    print("Recycle Length = %s " %(value(instance.LRe)))
for i in sorted(instance.Comp):
    print("Kap[%s] = %s" %(i,value(instance.Kap[i])))
    print("H[%s] = %s" %(i,value(instance.H[i])))
    if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
        print("b[%s] = %s" %(i,value(instance.b[i])))
    if Axial_D == True:
        print(f"Dax[{i}] = {value(instance.Dax[i])}")
print("--------------------------------------------------")
if Tikhonov == True:
    print("--------------------------------------------------")
    if EVM == True:
        print(f"regularization parameters = {instance.rho[1],instance.rho[2]}")
    elif EVM == False:
        print(f"regularization parameters = {instance.rho[1]}")
    print("--------------------------------------------------")
    if DV == True:
        print(f"Recycle Length Error = {value(instance.LRe)-L2LRe}")
    for i in sorted(instance.Comp):
        print(f"Kap[{i}] Error = {value(instance.Kap[i])-L2Kap[i-1]}")
        print(f"H[{i}] Error = {value(instance.H[i])-L2H[i-1]}")
        if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
            print(f"b[{i}] Error = {value(instance.b[i])-L2b[i-1]}")
        if Axial_D == True:
            print(f"Dax[{i}] Error = {value(instance.Dax[i])-L2Dax[i-1]}")
print("--------------------------------------------------")
if Tikhonov == True:
    if DV == True:
        print(f"Recycle Length Error [%] = {(value(instance.LRe)-L2LRe)*100/(L2LRe)}")
    for i in sorted(instance.Comp):
        print(f"Kap[{i}] Error [%] = {(value(instance.Kap[i])-L2Kap[i-1])*100/L2Kap[i-1]}")
        print(f"H[{i}] Error [%] = {(value(instance.H[i])-L2H[i-1])*100/L2H[i-1]}")
        if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
            print(f"b[{i}] Error [%] = {(value(instance.b[i])-L2b[i-1])*100/L2b[i-1]}")
        if Axial_D == True:
            print(f"Dax[{i}] Error = {(value(instance.Dax[i])-L2Dax[i-1])*100/L2Dax[i-1]}")
print("--------------------------------------------------")
if Tikhonov == True:
    print("L2 parameters in Tikhonov regularization term")
    if DV == True:
        print(f"Recycle Length = {L2LRe}")
    for i in sorted(instance.Comp):
        print("Kap[%s] = %s" %(i,L2Kap[i-1]))
        print("H[%s] = %s" %(i,L2H[i-1]))
        if Isotype == "Langmuir" or Isotype == "anti-Langmuir":
            print("b[%s] = %s" %(i,L2b[i-1]))
        if Axial_D == True:
            print(f"Dax[{i}] = {L2Dax[i-1]}")
print("--------------------------------------------------")

print("\n")

FittingTerm = np.zeros(5)
for d in instance.Data:
    FittingTerm[0] += value((instance.CbarE[d,1] - instance.averageCE[d,1])**2 + (instance.CbarR[d,1] - instance.averageCR[d,1])**2 + (instance.CbarE[d,2] - instance.averageCE[d,2])**2 + (instance.CbarR[d,2] - instance.averageCR[d,2])**2)
print("--------------------------------------------------")
print(f"Fitting term = {FittingTerm[0]}")
print(f"Fitting term divided by sample size = {(1/len(instance.Data))*(FittingTerm[0])}")
for d in instance.Data:
    FittingTerm[1] += value((instance.CbarE[d,1] - instance.averageCE[d,1])**2)
    FittingTerm[2] += value((instance.CbarR[d,1] - instance.averageCR[d,1])**2)
    FittingTerm[3] += value((instance.CbarE[d,2] - instance.averageCE[d,2])**2 )
    FittingTerm[4] += value((instance.CbarR[d,2] - instance.averageCR[d,2])**2)
    
print("--------------------------------------------------")
print(f"Fitting term of A in Extract = {FittingTerm[1]}")
print(f"Fitting term of A in Raffinate = {FittingTerm[2]}")
print(f"Fitting term of B in Extract = {FittingTerm[3]}")
print(f"Fitting term of B in Raffinate = {FittingTerm[4]}")
print("--------------------------------------------------")
print(f"Fitting term of A in Extract divided by sample size = {(1/len(instance.Data))*FittingTerm[1]}")
print(f"Fitting term of A in Raffinate divided by sample size = {(1/len(instance.Data))*FittingTerm[2]}")
print(f"Fitting term of B in Extract divided by sample size = {(1/len(instance.Data))*FittingTerm[3]}")
print(f"Fitting term of B in Raffinate divided by sample size = {(1/len(instance.Data))*FittingTerm[4]}")
print("--------------------------------------------------")

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
    if ini_plots == True:
        plt.figure(i+1+NData)
    elif ini_plots == False:
        plt.figure(i+1)
    else:
        print("ERROR: Confirm ini_plots definition")
    plt.plot(X, CA[i], 'o-', label='Comp A')
    plt.plot(X, CB[i], 'x-', label='Comp B')
    plt.legend()
    plt.xlim(0, ColL*Colum)
    plt.ylim(0, 1.2)
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
        plt.ylim(0, 1.2)
        plt.xlim(0, ColL*Colum)
        plt.plot(x,u1[i,:], 'o-', label='Comp A')
        plt.plot(x,u2[i,:], 'x-', label='Comp B')
        plt.xlabel('x [m]')
        plt.ylabel('C/C(feed) [-]')
        plt.title(f"Data{l+1} SMB Concentration profile")
        plt.legend()


    nframe = Nfet*value(instance.NCP)+1
    nframe = nframe*4
    ani = animation.FuncAnimation(fig,update,fargs = (C1_ani,C2_ani,X,'Wave motion'),interval=100,frames=nframe)
    filename = dirname02 + f"Data{l+1}_Concentration_Animation.gif"
    ani.save(filename, writer="pillow", dpi=500)