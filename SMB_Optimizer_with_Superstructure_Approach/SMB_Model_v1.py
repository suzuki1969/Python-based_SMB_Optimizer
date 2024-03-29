# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kosuke Sekishita and Kensuke Suzuki in 2021
# The method used in this file is the same as study by Kawajiri and Biegler 
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.
# The superstructuer approach used in this script comes from Sreedhar and Kawajiri
# Bibliographic information: Sreedhar B, Kawajiri Y. Multi-column chromatographic process development using simulated moving bed superstructure and simultaneous optimization - model correction framework. Chem Eng Sci. 2014;116:428-441

# ==================================================================
# Generalized SMB optimization
# ==================================================================

from pyomo.environ import *
from pyomo.dae import *
from Initdata_SMB_v1 import *
import numpy as np

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
m.NFEXX = Param(initialize = Nfex-1)
m.NCP = Param(initialize = NCP)
m.NComp = Param(initialize = NComp)

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

def FEXX_init(m):
	retval = []
	j = 1
	while(j<=m.NFEXX):
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

# -------------------------------------------------------------------
# Indices / Dimensions
# -------------------------------------------------------------------

m.FET = Set(initialize = FET_init, ordered = True)
m.Col = Set(initialize = Col_init, ordered = True, dimen=1)
m.FEX = Set(initialize = FEX_init, ordered = True)
m.FEXX = Set(initialize = FEXX_init, ordered = True)
m.CP = Set(initialize = CP_init, ordered = True)
m.Comp = Set(initialize = Comp_init, ordered = True, dimen=1)

# -------------------------------------------------------------------
# Physical Paramters
# -------------------------------------------------------------------

m.L = Param(initialize = ColL)
m.Lmax = Param(initialize = m.L*1.2)
m.Lmin = Param(initialize = m.L/1.2)

def L_init(m,i):
	return m.L

m.Lvar = Var(m.Col, initialize = L_init, bounds = (m.Lmin, m.Lmax))

m.eb = Param(initialize = eb)
m.CF = Param(m.Comp, initialize = CFeed) # feed concentration

m.K = Param(m.Comp, initialize = K) # Henry's constant
m.Kap = Param(m.Comp, initialize = Kap) # overall mass transfer coefficient

if Isotype == "Henry":
	pass
elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
	m.b = Param(m.Comp, initialize = b) # affinity parameter
else:
	print("ERROR: Confirm Isotype definition")

# -------------------------------------------------------------------
# Product Specification Parameters
# -------------------------------------------------------------------

m.ExtractRecMin = Param(m.Comp, initialize=ExtractRecMin)
m.ExtractPurMin = Param(m.Comp, initialize=ExtractPurMin)
m.RaffinateRecMin = Param(m.Comp, initialize=RaffinateRecMin)
m.RaffinatePurMin = Param(m.Comp, initialize=RaffinatePurMin)
m.ConcentrationMin = Param(m.Comp, initialize=ConcentrationMin)

# -------------------------------------------------------------------
# State Variables
# -------------------------------------------------------------------

m.x = ContinuousSet(bounds=(0,m.L))

m.t = ContinuousSet(bounds=(0,1))

m.HT = Var(m.t, bounds = ((0.5/m.NFET), (2/m.NFET)))

m.C = Var(m.Comp, m.Col, m.t, m.x) # concentration in liquid phase
m.Q = Var(m.Comp, m.Col, m.t, m.x) # concentration in solid phase
m.CEQ = Var(m.Comp, m.Col, m.t, m.x) # equilibriun concentration in liquid phase

m.intCR = Var(m.Comp, m.Col)
m.intCE = Var(m.Comp, m.Col)

if DV == 'yes':
	m.xRe = ContinuousSet(bounds=(0,1)) # m.x in dead volume
	m.CRe = Var(m.Comp, m.t, m.xRe) # m.C in dead volume
elif DV == 'no':
    pass
else:
    print("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# Control Variables
# -------------------------------------------------------------------

m.U = Var(m.Col, m.t, within=PositiveReals, initialize=1.0) # fluid velocity in column
m.UF = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Feed
m.UD = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Desorbent
m.UR = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Raffinate
m.UE = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Extract

if DV == 'yes':
	m.URe = Var(m.t, within=PositiveReals) # fluid velocity in dead volume
	m.LRe = Param(initialize = DeadVolume) # length of dead volume
elif DV == 'no':
    pass
else:
    print("ERROR: Confirm DeadVolume definition")

m.StepTime = Var(bounds = (5,3600))

# -------------------------------------------------------------------
# State Derivatives
# -------------------------------------------------------------------

m.dCdt = DerivativeVar(m.C, wrt=m.t)
m.dQdt = DerivativeVar(m.Q, wrt=m.t)
m.dCdx = DerivativeVar(m.C, wrt=m.x)

if DV == 'yes':
	m.dCRedt = DerivativeVar(m.CRe, wrt=m.t)
	m.dCRedx = DerivativeVar(m.CRe, wrt=m.xRe)
elif DV == 'no':
    pass
else:
    print("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# Flow rate balances
# -------------------------------------------------------------------

def FlowBalance_rule(m,j,k):
	if j==m.NColumn: return Constraint.Skip
	else: return m.U[j+1,k]==m.U[j,k]-m.UR[j,k]-m.UE[j,k]+m.UD[j+1,k]+m.UF[j+1,k]

m.FlowBalance=Constraint(m.Col, m.t, rule=FlowBalance_rule)

if DV == 'yes':
	def FlowBalanceD1_rule(m,k):
		return m.U[1,k] == m.URe[k] + m.UD[1,k] + m.UF[1,k]
	m.FlowBalanceD1=Constraint(m.t, rule=FlowBalanceD1_rule)

	# def FlowBalance4D_rule(m,k):
	# 	return m.URe[k] == m.U[Colum,k] - m.UR[Colum,k] - m.UE[Colum,k]
	# m.FlowBalance4D=Constraint(m.t, rule=FlowBalance4D_rule)
elif DV == 'no':
	pass
else:
    print("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# Boundary Conditions
# -------------------------------------------------------------------

def BoundaryCondition_rule(m, i,j,k):
   if k == 0.0: return Constraint.Skip
   if j==m.NColumn: return Constraint.Skip
   else: return m.C[i,j+1,k,0]*m.U[j+1,k]==m.C[i,j,k,m.L]*(m.U[j,k] - m.UE[j,k] - m.UR[j,k]) + m.CF[i]*m.UF[j+1,k]

m.BoundaryCondition=Constraint(m.Comp, m.Col, m.t, rule=BoundaryCondition_rule)

if DV == 'yes':
	def BoundaryCondition4D_rule(m,i,k):
		if k == 0 : return Constraint.Skip
		return m.CRe[i,k,0]*m.URe[k]==m.C[i,Colum,k,m.L]*(m.U[Colum,k] - m.UE[Colum,k] - m.UR[Colum,k])
	m.BoundaryCondition4D=Constraint(m.Comp, m.t, rule=BoundaryCondition4D_rule)

	def BoundaryConditionD1_rule(m,i,k):
		if k == 0 : return Constraint.Skip
		return m.C[i,1,k,0]*m.U[1,k]==m.CRe[i,k,1]*m.URe[k] + m.CF[i]*m.UF[1,k]
	m.BoundaryConditionD1=Constraint(m.Comp, m.t, rule=BoundaryConditionD1_rule)
elif DV == 'no':
	def BoundaryCondition1_rule(m,i,k):
		if k == 0.0: return Constraint.Skip
		return m.C[i,1,k,0]*m.U[1,k]==m.C[i,Colum,k,m.L]*(m.U[Colum,k]-m.UE[Colum,k]-m.UR[Colum,k]) + m.CF[i]*m.UF[1,k]
	m.BoundaryCondition1=Constraint(m.Comp, m.t, rule=BoundaryCondition1_rule)
else:
    print("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# Model equations
# -------------------------------------------------------------------

def MassBalanceLiquid_rule(m, i,j,k,l):
	if l==0: return Constraint.Skip
	else: return m.eb*m.dCdt[i,j,k,l]/(m.HT[k]*m.NFET)/m.StepTime + (1-m.eb)*m.dQdt[i,j,k,l]/(m.HT[k]*m.NFET)/m.StepTime + m.U[j,k]*m.dCdx[i,j,k,l] == 0

def MassBalanceSolid_rule(m, i,j,k,l):
	return (1-m.eb)*m.dQdt[i,j,k,l]/(m.HT[k]*m.NFET)/m.StepTime == m.Kap[i]*(m.C[i,j,k,l]-m.CEQ[i,j,k,l])

if Isotype == "Henry":
	def Equilibrium_rule(m, i,j,k,l):
		return m.Q[i,j,k,l] == m.K[i]*m.CEQ[i,j,k,l]
elif Isotype == "Langmuir":
	def Equilibrium_rule(m, i,j,k,l):
		return m.Q[i,j,k,l]*(1.0 + m.b[1]*m.CEQ[1,j,k,l] + m.b[2]*m.CEQ[2,j,k,l]) == m.K[i]*m.CEQ[i,j,k,l]
elif Isotype == "anti-Langmuir":
	def Equilibrium_rule(m, i,j,k,l):
		return m.Q[i,j,k,l]*(1.0 - m.b[1]*m.CEQ[1,j,k,l] - m.b[2]*m.CEQ[2,j,k,l]) == m.K[i]*m.CEQ[i,j,k,l]
else:
	print("ERROR: Confirm Isotype definition")

m.MassBalanceLiquid=Constraint(m.Comp, m.Col, m.t, m.x, rule=MassBalanceLiquid_rule)
m.MassBalanceSolid=Constraint(m.Comp, m.Col, m.t, m.x, rule=MassBalanceSolid_rule)
m.Equilibrium=Constraint(m.Comp, m.Col, m.t, m.x, rule=Equilibrium_rule)

if DV == 'yes':
	def DeadVolume_rule(m, i,j,k):
		if k==0: return Constraint.Skip
		else: return m.dCRedt[i,j,k]/(m.StepTime*m.HT[j]*Nfet) + m.URe[j]*m.dCRedx[i,j,k]/m.LRe == 0
	m.DeadVolume=Constraint(m.Comp, m.t, m.xRe, rule=DeadVolume_rule)
elif DV == 'no':
    pass
else:
    print("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# No Backflow constriant
# -------------------------------------------------------------------

def NoBackFlow_rule(m, j, k):
	return m.U[j,k]-m.UD[j,k]-m.UF[j,k] >=0

m.NoBackFlow=Constraint(m.Col, m.t, rule=NoBackFlow_rule)

# -------------------------------------------------------------------
# Cyclic Steady State Convergence
# -------------------------------------------------------------------

def C_CSS_p_rule(m, i,j,p):
	if j==1: return Constraint.Skip
	else: return m.C[i,j,1,p] - m.C[i,j-1,0,p] == 0

def Q_CSS_p_rule(m, i,j,p):
	if j==1: return Constraint.Skip
	else: return m.Q[i,j,1,p] - m.Q[i,j-1,0,p] == 0

def C_CSS_1_p_rule(m, i,p):
	return m.C[i,1,1,p] - m.C[i,Colum,0,p]==0

def Q_CSS_1_p_rule(m, i,p):
	return m.Q[i,1,1,p] - m.Q[i,Colum,0,p]==0

m.C_CSS_p=Constraint(m.Comp, m.Col, m.x, rule=C_CSS_p_rule)
m.Q_CSS_p=Constraint(m.Comp, m.Col, m.x, rule=Q_CSS_p_rule)
m.C_CSS_1_p=Constraint(m.Comp, m.x, rule=C_CSS_1_p_rule)
m.Q_CSS_1_p=Constraint(m.Comp, m.x, rule=Q_CSS_1_p_rule)

if DV == 'yes':
	def CRe_CSS_p_rule(m, i,p):
		return m.CRe[i,1,p] - m.CRe[i,0,p]==0
	m.CRe_CSS_p=Constraint(m.Comp, m.xRe, rule=CRe_CSS_p_rule)
elif DV == 'no':
    pass
else:
    print("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# Purity and Recovery Constraints
# -------------------------------------------------------------------

def ExtractRecoveryConstraint_rule(m):
    return sum(m.intCE[2,j] for j in m.Col) >= (m.ExtractRecMin[2]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[2]/m.NCP for j in m.Col for k in m.t)

def ExtractPurityConstraint_rule(m):
	return sum(m.intCE[2,j] for j in m.Col) >= (m.ExtractPurMin[2]*((dS+100)/100)/100)*(sum(m.intCE[n,j] for n in m.Comp for j in m.Col))

def RaffinateRecoveryConstraint_rule(m):
	return sum(m.intCR[1,j] for j in m.Col) >= (m.RaffinateRecMin[1]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[1]/m.NCP for j in m.Col for k in m.t)

def RaffinatePurityConstraint_rule(m):
	return sum(m.intCR[1,j] for j in m.Col) >= (m.RaffinatePurMin[1]*((dS+100)/100)/100)*(sum(m.intCR[n,j] for n in m.Comp for j in m.Col))

m.ExtractRecoveryConstraint = Constraint(rule=ExtractRecoveryConstraint_rule)
m.ExtractPurityConstraint = Constraint(rule=ExtractPurityConstraint_rule)
m.RaffinateRecoveryConstraint = Constraint(rule=RaffinateRecoveryConstraint_rule)
m.RaffinatePurityConstraint = Constraint(rule=RaffinatePurityConstraint_rule)

"""
# -------------------------------------------------------------------
# Purity
# -------------------------------------------------------------------
m.p_EB = Var()
@m.Constraint()
def purity_EB(m):
    return sum(m.intCE[2,j] for j in m.Col) == (m.p_EB/100)*(sum(m.intCE[n,j] for n in m.Comp for j in m.Col))

m.p_RA = Var()
@m.Constraint()
def purity_RA(m):
    return sum(m.intCE[1,j] for j in m.Col) == (m.p_RA/100)*(sum(m.intCE[n,j] for n in m.Comp for j in m.Col))
"""