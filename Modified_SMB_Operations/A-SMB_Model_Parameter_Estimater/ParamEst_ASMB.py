# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# ASMB model was constructed by Hideki Harada and Kensuke Suzuki in 2020
# The mathematical model performed in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

# ==================================================================
# Generalized SMB optimization
# ==================================================================

from pyomo.environ import *
from pyomo.dae import ContinuousSet, DerivativeVar
from ParamEst_Initdata_ASMB import *

Nfex=Nfex+1

m = AbstractModel()

# -------------------------------------------------------------------
# PARAMETERS 
# NFET = Number of Finite Elements in Time
# NColumn = Number of Columns
# NFEX = Number of Finite Elements in Space
# NCP = Number of Collocation Points
# NComp = Number of Components 
# -------------------------------------------------------------------

m.NFET = Param(initialize = Nfet) 
m.NColumn = Param(initialize = Colum) 
m.NFEX = Param(initialize = Nfex)     
m.NCP = Param(initialize = NCP)
m.NComp = Param(initialize = NComp)

m.NData = Param(initialize = NData)


def FET_init(m):
	retval = []
	j = 1
	while(j<=m.NFET):
		retval.append(j)
		j=j+1
	return retval

def Col_init(m):
	retval = []
	j = 1
	while(j<=m.NColumn):
		retval.append(j)
		j=j+1
	return retval

def FEX_init(m):
	retval = []
	j = 1
	while(j<=m.NFEX):
		retval.append(j)
		j=j+1
	return retval

def CP_init(m):
	retval = []
	j = 1
	while(j<=m.NCP):
		retval.append(j)
		j=j+1
	return retval

def Comp_init(m):
	retval = []
	j = 1
	while(j<=m.NComp):
		retval.append(j)
		j=j+1
	return retval

def Data_init(m):
	retval = []
	j = 1
	while(j<=m.NData):
		retval.append(j)
		j=j+1
	return retval

# -------------------------------------------------------------------
# Indices / Dimensions 
# -------------------------------------------------------------------

m.FET = Set(initialize = FET_init, ordered = True)  
m.Col = Set(initialize = Col_init, ordered = True, dimen=1)  
m.FEX = Set(initialize = FEX_init, ordered = True)  
m.CP = Set(initialize = CP_init, ordered = True)  
m.Comp = Set(initialize = Comp_init, ordered = True, dimen=1)  
m.Ridge = Set(initialize = [1,2])  
m.Data = Set(initialize = Data_init, ordered = True, dimen=1)

# -------------------------------------------------------------------
# Physical Paramters
# -------------------------------------------------------------------

m.L = Param(initialize = ColL) # column length

m.eb = Param(initialize = eb) # void fraction
m.CF = Param(m.Data, m.Comp, initialize = CFeed) # feed concentration
m.K = Var(m.Comp, within=PositiveReals) # Henry's constant
m.Kap = Var(m.Comp, within=PositiveReals) # overall mass transfer coefficient

if Isotype == "Henry":
	pass
elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
	m.b = Var(m.Comp, within=PositiveReals) # affinity parameter
else:
	print("ERROR: Confirm Isotype definition")

m.averageCE = Var(m.Data, m.Comp) # average conc in Extract (model)
m.averageCR = Var(m.Data, m.Comp) # average conc in Raffinate (model)

m.CbarE = Param(m.Data, m.Comp, initialize = CExtract) # average conc in Extract (experiment)
m.CbarR = Param(m.Data, m.Comp, initialize = CRaffinate) # average conc in Raffinate (experiment)

m.rho = Param(m.Ridge, initialize = Rho) # regularization parameters

data = DataPortal(model=m)

# -------------------------------------------------------------------
# State Variables
# -------------------------------------------------------------------

m.t = ContinuousSet(bounds=(0,1)) # time

m.HT = Var(m.Data, m.t, bounds = ((0.5/m.NFET), (2/m.NFET))) # length of a finite element in time

m.x = ContinuousSet(bounds=(0,m.L)) # space

m.C = Var(m.Data, m.Comp, m.Col, m.t, m.x) # concentration in liquid phase
m.Q = Var(m.Data, m.Comp, m.Col, m.t, m.x) # concentration in solid phase
m.CEQ = Var(m.Data, m.Comp, m.Col, m.t, m.x) # equilibriun concentration in liquid phase

m.xRe = ContinuousSet(bounds=(0,1)) # m.x in dead volume
m.CRe = Var(m.Data, m.Comp, m.t, m.xRe) # m.C in dead volume

# -------------------------------------------------------------------
# Control Variables
# -------------------------------------------------------------------

m.U = Var(m.Data, m.Col, m.t, within=PositiveReals ) # fluid velocity in column
m.UF = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Feed
m.UD = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Desorbent
m.UR = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Raffinate
m.UE = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Extract

m.URe = Var(m.Data, m.t, within=PositiveReals) # fluid velocity in dead volume
m.LRe = Var(initialize = 0.05*Colum*ColL, within=PositiveReals) # length of dead volume

m.StepTime = Var(m.Data, bounds = (5,3600))

# -------------------------------------------------------------------
# State Derivatives
# -------------------------------------------------------------------

m.dCdt = DerivativeVar(m.C, wrt=m.t)
m.dQdt = DerivativeVar(m.Q, wrt=m.t)
m.dCdx = DerivativeVar(m.C, wrt=m.x)

m.dCRedt = DerivativeVar(m.CRe, wrt=m.t)
m.dCRedx = DerivativeVar(m.CRe, wrt=m.xRe)

# -------------------------------------------------------------------
# Flow rate balances
# -------------------------------------------------------------------

def FlowBalance_rule(m,d,j,k):
	if j==m.NColumn: return Constraint.Skip
	else: return m.U[d,j+1,k]==m.U[d,j,k]-m.UR[d,j,k]-m.UE[d,j,k]+m.UD[d,j+1,k]+m.UF[d,j+1,k]

m.FlowBalance=Constraint(m.Data, m.Col, m.t, rule=FlowBalance_rule)

def FlowBalance4D_rule(m,d,k):
    return m.URe[d,k] == m.U[d,Colum,k]
m.FlowBalance4D=Constraint(m.Data, m.t, rule=FlowBalance4D_rule)

def FlowBalanceD1_rule(m,d,k):
    return m.U[d,1,k] == m.URe[d,k] + m.UD[d,1,k]
m.FlowBalanceD1=Constraint(m.Data, m.t, rule=FlowBalanceD1_rule)

# -------------------------------------------------------------------
# Boundary Conditions
# -------------------------------------------------------------------

def BoundaryCondition_rule(m, d,i,j,k):
    if k == 0 : return Constraint.Skip
    if j==m.NColumn: return Constraint.Skip
    else: return m.C[d,i,j+1,k,0]*m.U[d,j+1,k]==m.C[d,i,j,k,m.L]*(m.U[d,j,k] - m.UE[d,j,k] - m.UR[d,j,k]) + m.CF[d,i]*m.UF[d,j+1,k]
m.BoundaryCondition=Constraint(m.Data, m.Comp, m.Col, m.t, rule=BoundaryCondition_rule)

def BoundaryCondition4D_rule(m,d,i,k):
    if k == 0 : return Constraint.Skip 
    return m.CRe[d,i,k,0]*m.URe[d,k]==m.C[d,i,Colum,k,m.L]*m.U[d,Colum,k] 
m.BoundaryCondition4D=Constraint(m.Data, m.Comp, m.t, rule=BoundaryCondition4D_rule)

def BoundaryConditionD1_rule(m,d,i,k):
    if k == 0 : return Constraint.Skip 
    return m.C[d,i,1,k,0]*m.U[d,1,k]==m.CRe[d,i,k,1]*m.URe[d,k] 
m.BoundaryConditionD1=Constraint(m.Data, m.Comp, m.t, rule=BoundaryConditionD1_rule)

# -------------------------------------------------------------------
# Model equations
# -------------------------------------------------------------------

def MassBalanceLiquid_rule(m, d,i,j,k,l):
	if l==0: return Constraint.Skip
	else: return m.eb*m.dCdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*Nfet) + (1-m.eb)*m.dQdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*Nfet) + m.U[d,j,k]*m.dCdx[d,i,j,k,l] == 0 

def MassBalanceSolid_rule(m, d,i,j,k,l):
    return (1.0-m.eb)*m.dQdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*Nfet) == m.Kap[i]*(m.C[d,i,j,k,l]-m.CEQ[d,i,j,k,l])

if Isotype == "Henry":
	def Equilibrium_rule(m, d,i,j,k,l):
		return m.Q[d,i,j,k,l] == m.K[i]*m.CEQ[d,i,j,k,l]
elif Isotype == "Langmuir":
	def Equilibrium_rule(m, d,i,j,k,l):
		return m.Q[d,i,j,k,l]*(1.0 + m.b[1]*m.CEQ[d,1,j,k,l] + m.b[2]*m.CEQ[d,2,j,k,l]) == m.K[i]*m.CEQ[d,i,j,k,l]
elif Isotype == "anti-Langmuir":
	def Equilibrium_rule(m, d,i,j,k,l):
		return m.Q[d,i,j,k,l]*(1.0 - m.b[1]*m.CEQ[d,1,j,k,l] - m.b[2]*m.CEQ[d,2,j,k,l]) == m.K[i]*m.CEQ[d,i,j,k,l]
else:
	print("ERROR: Confirm Isotype definition")

m.MassBalanceLiquid=Constraint(m.Data, m.Comp, m.Col, m.t, m.x, rule=MassBalanceLiquid_rule)
m.MassBalanceSolid=Constraint(m.Data, m.Comp, m.Col, m.t, m.x, rule=MassBalanceSolid_rule)
m.Equilibrium=Constraint(m.Data, m.Comp, m.Col, m.t, m.x, rule=Equilibrium_rule)

def DeadVolume_rule(m, d,i,j,k):
    if k==0: return Constraint.Skip
    else: return m.dCRedt[d,i,j,k]/(m.StepTime[d]*m.HT[d,j]*Nfet) + m.URe[d,j]*m.dCRedx[d,i,j,k]/m.LRe == 0
m.DeadVolume=Constraint(m.Data, m.Comp, m.t, m.xRe, rule=DeadVolume_rule)

# -------------------------------------------------------------------
# No Backflow constriant
# -------------------------------------------------------------------

def NoBackFlow_rule(m, d,j,k):
	return m.U[d,j,k]-m.UD[d,j,k]-m.UF[d,j,k] >=0
m.NoBackFlow=Constraint(m.Data, m.Col, m.t, rule=NoBackFlow_rule)

# -------------------------------------------------------------------
# Cyclic Steady State Convergence
# -------------------------------------------------------------------

def C_CSS_p_rule(m, d,i,j,p):
	if j==1: return Constraint.Skip
	else: return m.C[d,i,j,1,p] - m.C[d,i,j-1,0,p] == 0

def Q_CSS_p_rule(m, d,i,j,p):
	if j==1: return Constraint.Skip
	else: return m.Q[d,i,j,1,p] - m.Q[d,i,j-1,0,p] == 0

def C_CSS_1_p_rule(m, d,i,p):
	return m.C[d,i,1,1,p] - m.C[d,i,Colum,0,p]==0

def Q_CSS_1_p_rule(m, d,i,p):
	return m.Q[d,i,1,1,p] - m.Q[d,i,Colum,0,p]==0

m.C_CSS_p=Constraint(m.Data, m.Comp, m.Col, m.x, rule=C_CSS_p_rule)
m.Q_CSS_p=Constraint(m.Data, m.Comp, m.Col, m.x, rule=Q_CSS_p_rule)
m.C_CSS_1_p=Constraint(m.Data, m.Comp, m.x, rule=C_CSS_1_p_rule)
m.Q_CSS_1_p=Constraint(m.Data, m.Comp, m.x, rule=Q_CSS_1_p_rule)

def CRe_CSS_p_rule(m, d,i,p):
	return m.CRe[d,i,1,p] - m.CRe[d,i,0,p]==0
m.CRe_CSS_p=Constraint(m.Data, m.Comp, m.xRe, rule=CRe_CSS_p_rule)
