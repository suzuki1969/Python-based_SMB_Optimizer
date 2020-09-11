# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kensuke Suzuki in 2020
# The mathematical model performed in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.

# ==================================================================
# Generalized SMB optimization
# ==================================================================

from pyomo.environ import *
from pyomo.dae import *
from Initdata_SMB_Optimization import *
import numpy as np
import sys

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
m.Col = Set(initialize = Col_init, ordered = True)
m.FEX = Set(initialize = FEX_init, ordered = True)
m.FEXX = Set(initialize = FEXX_init, ordered = True)
m.CP = Set(initialize = CP_init, ordered = True)
m.Comp = Set(initialize = Comp_init, ordered = True)

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

m.H = Param(m.Comp, initialize = H) # Henry's constant
m.Kap = Param(m.Comp, initialize = Kap) # overall mass transfer coefficient

if Isotype == "Henry":
	pass
elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
	m.b = Param(m.Comp, initialize = b) # affinity parameter
else:
	sys.exit("ERROR: Confirm Isotype definition")

if Axial_D == True:
	m.Dax = Param(m.Comp, initialize = Dax) # Axial Dispersion Coefficient
elif Axial_D == False:
	pass
else:
	sys.exit("ERROR: Confirm Axial_D definition")

# -------------------------------------------------------------------
# Product Specification Parameters
# -------------------------------------------------------------------

m.ExtractRecMin = Param(m.Comp, initialize=ExtractRecMin)
m.ExtractPurMin = Param(m.Comp, initialize=ExtractPurMin)
m.RaffinateRecMin = Param(m.Comp, initialize=RaffinateRecMin)
m.RaffinatePurMin = Param(m.Comp, initialize=RaffinatePurMin)

# -------------------------------------------------------------------
# State Variables
# -------------------------------------------------------------------

m.x = ContinuousSet(bounds=(0,m.L))

m.t = ContinuousSet(bounds=(0,1))

m.HT = Var(m.t, bounds = ((0.5/m.NFET), (2/m.NFET)))

m.C = Var(m.Comp, m.Col, m.t, m.x) # concentration in liquid phase
m.Q = Var(m.Comp, m.Col, m.t, m.x) # concentration in solid phase

m.intCR = Var(m.Comp, m.Col)
m.intCE = Var(m.Comp, m.Col)

if Phase == 'Liquid':
	m.CEQ = Var(m.Comp, m.Col, m.t, m.x) # equilibriun concentration in liquid phase
elif Phase == 'Solid':
	m.QEQ = Var(m.Comp, m.Col, m.t, m.x) # equilibriun concentration in sloid phase
else:
	sys.exit("ERROR: Confirm Phase definition")

if DV == True:
	m.xRe = ContinuousSet(bounds=(0,1)) # m.x in dead volume
	m.CRe = Var(m.Comp, m.t, m.xRe) # m.C in dead volume
elif DV == False:
    pass
else:
    sys.exit("ERROR: Confirm Dead Volume definition")

# -------------------------------------------------------------------
# Control Variables
# -------------------------------------------------------------------

m.U = Var(m.Col, m.t, within=PositiveReals, initialize=1.0) # fluid velocity in column
m.UF = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Feed
m.UD = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Desorbent
m.UR = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Raffinate
m.UE = Var(m.Col, m.t, within=PositiveReals) # fluid velocity of Extract

if DV == True:
	m.URe = Var(m.t, within=PositiveReals) # fluid velocity in dead volume
	m.LRe = Param(initialize = DeadVolume) # length of dead volume

m.StepTime = Var(bounds = (5,3600))

# -------------------------------------------------------------------
# State Derivatives
# -------------------------------------------------------------------

m.dCdt = DerivativeVar(m.C, wrt=m.t)
m.dQdt = DerivativeVar(m.Q, wrt=m.t)
m.dCdx = DerivativeVar(m.C, wrt=m.x)

if Axial_D == True:
	m.dCdx2 = DerivativeVar(m.C, wrt=(m.x, m.x)) # second-order differentiated liquid concentration with respect to space

if DV == True:
	m.dCRedt = DerivativeVar(m.CRe, wrt=m.t)
	m.dCRedx = DerivativeVar(m.CRe, wrt=m.xRe)

# -------------------------------------------------------------------
# Flow rate balances
# -------------------------------------------------------------------

def FlowBalance_rule(m,j,k):
	if j==m.NColumn: return Constraint.Skip
	else: return m.U[j+1,k]==m.U[j,k]-m.UR[j,k]-m.UE[j,k]+m.UD[j+1,k]+m.UF[j+1,k]

m.FlowBalance=Constraint(m.Col, m.t, rule=FlowBalance_rule)

if DV == True:
	def FlowBalanceD1_rule(m,k):
		return m.U[1,k] == m.URe[k] + m.UD[1,k] + m.UF[1,k]
	m.FlowBalanceD1=Constraint(m.t, rule=FlowBalanceD1_rule)

# -------------------------------------------------------------------
# Boundary Conditions
# -------------------------------------------------------------------

def BoundaryCondition_rule(m, i,j,k):
   if k == 0.0: return Constraint.Skip
   if j==m.NColumn: return Constraint.Skip
   else: return m.C[i,j+1,k,0]*m.U[j+1,k]==m.C[i,j,k,m.L]*(m.U[j,k] - m.UE[j,k] - m.UR[j,k]) + m.CF[i]*m.UF[j+1,k]

m.BoundaryCondition=Constraint(m.Comp, m.Col, m.t, rule=BoundaryCondition_rule)

if DV == True:
	def BoundaryCondition4D_rule(m,i,k):
		if k == 0 : return Constraint.Skip 
		return m.CRe[i,k,0]*m.URe[k]==m.C[i,Colum,k,m.L]*(m.U[Colum,k] - m.UE[Colum,k] - m.UR[Colum,k])
	m.BoundaryCondition4D=Constraint(m.Comp, m.t, rule=BoundaryCondition4D_rule)

	def BoundaryConditionD1_rule(m,i,k):
		if k == 0 : return Constraint.Skip 
		return m.C[i,1,k,0]*m.U[1,k]==m.CRe[i,k,1]*m.URe[k] + m.CF[i]*m.UF[1,k]
	m.BoundaryConditionD1=Constraint(m.Comp, m.t, rule=BoundaryConditionD1_rule)
elif DV == False:
	def BoundaryCondition1_rule(m,i,k):
		if k == 0.0: return Constraint.Skip
		return m.C[i,1,k,0]*m.U[1,k]==m.C[i,Colum,k,m.L]*(m.U[Colum,k]-m.UE[Colum,k]-m.UR[Colum,k]) + m.CF[i]*m.UF[1,k]
	m.BoundaryCondition1=Constraint(m.Comp, m.t, rule=BoundaryCondition1_rule)

# -------------------------------------------------------------------
# Model equations
# -------------------------------------------------------------------

if Axial_D == True:
	def MassBalanceLiquid_rule(m, i,j,k,l):
		if l==0: return Constraint.Skip
		else: return m.eb*m.dCdt[i,j,k,l]/(m.StepTime*m.HT[k]*m.NFET) + (1-m.eb)*m.dQdt[i,j,k,l]/(m.StepTime*m.HT[k]*m.NFET) + m.U[j,k]*m.dCdx[i,j,k,l] == m.Dax[i]*m.dCdx2[i,j,k,l] 
elif Axial_D == False:
	def MassBalanceLiquid_rule(m, i,j,k,l):
		if l==0: return Constraint.Skip
		else: return m.eb*m.dCdt[i,j,k,l]/(m.StepTime*m.HT[k]*m.NFET) + (1-m.eb)*m.dQdt[i,j,k,l]/(m.StepTime*m.HT[k]*m.NFET) + m.U[j,k]*m.dCdx[i,j,k,l] == 0 

if Phase == 'Liquid':
	def MassBalanceSolid_rule(m, i,j,k,l):
		return (1.0-m.eb)*m.dQdt[i,j,k,l]/(m.StepTime*m.HT[k]*m.NFET) == m.Kap[i]*(m.C[i,j,k,l]-m.CEQ[i,j,k,l])

	if Isotype == "Henry":
		def Equilibrium_rule(m, i,j,k,l):
			return m.Q[i,j,k,l] == m.H[i]*m.CEQ[i,j,k,l]
	elif Isotype == "Langmuir":
		def Equilibrium_rule(m, i,j,k,l):
			return m.Q[i,j,k,l]*(1.0 + m.b[1]*m.CEQ[1,j,k,l] + m.b[2]*m.CEQ[2,j,k,l]) == m.H[i]*m.CEQ[i,j,k,l]
	elif Isotype == "anti-Langmuir":
		def Equilibrium_rule(m, i,j,k,l):
			return m.Q[i,j,k,l]*(1.0 - m.b[1]*m.CEQ[1,j,k,l] - m.b[2]*m.CEQ[2,j,k,l]) == m.H[i]*m.CEQ[i,j,k,l]

elif Phase == 'Solid':
	def MassBalanceSolid_rule(m, i,j,k,l):
		return m.dQdt[i,j,k,l]/(m.StepTime*m.HT[k]*m.NFET) == m.Kap[i]*(m.QEQ[i,j,k,l]-m.Q[i,j,k,l])

	if Isotype == "Henry":
		def Equilibrium_rule(m, i,j,k,l):
			return m.QEQ[i,j,k,l] == m.H[i]*m.C[i,j,k,l]
	elif Isotype == "Langmuir":
		def Equilibrium_rule(m, i,j,k,l):
			return m.QEQ[i,j,k,l]*(1.0 + m.b[1]*m.C[1,j,k,l] + m.b[2]*m.C[2,j,k,l]) == m.H[i]*m.C[i,j,k,l]
	elif Isotype == "anti-Langmuir":
		def Equilibrium_rule(m, i,j,k,l):
			return m.QEQ[i,j,k,l]*(1.0 - m.b[1]*m.C[1,j,k,l] - m.b[2]*m.C[2,j,k,l]) == m.H[i]*m.C[i,j,k,l]

# def MassBalanceLiquid_rule(m, i,j,k,l):
# 	if l==0: return Constraint.Skip
# 	else: return m.eb*m.dCdt[i,j,k,l]/(m.HT[k]*m.NFET)/m.StepTime + (1-m.eb)*m.dQdt[i,j,k,l]/(m.HT[k]*m.NFET)/m.StepTime + m.U[j,k]*m.dCdx[i,j,k,l] == 0 

# def MassBalanceSolid_rule(m, i,j,k,l):
# 	return (1-m.eb)*m.dQdt[i,j,k,l]/(m.HT[k]*m.NFET)/m.StepTime == m.Kap[i]*(m.C[i,j,k,l]-m.CEQ[i,j,k,l]) 

# if Isotype == "Henry":
# 	def Equilibrium_rule(m, i,j,k,l):
# 		return m.Q[i,j,k,l] == m.K[i]*m.CEQ[i,j,k,l]
# elif Isotype == "Langmuir":
# 	def Equilibrium_rule(m, i,j,k,l):
# 		return m.Q[i,j,k,l]*(1.0 + m.b[1]*m.CEQ[1,j,k,l] + m.b[2]*m.CEQ[2,j,k,l]) == m.K[i]*m.CEQ[i,j,k,l]
# elif Isotype == "anti-Langmuir":
# 	def Equilibrium_rule(m, i,j,k,l):
# 		return m.Q[i,j,k,l]*(1.0 - m.b[1]*m.CEQ[1,j,k,l] - m.b[2]*m.CEQ[2,j,k,l]) == m.K[i]*m.CEQ[i,j,k,l]

m.MassBalanceLiquid=Constraint(m.Comp, m.Col, m.t, m.x, rule=MassBalanceLiquid_rule)
m.MassBalanceSolid=Constraint(m.Comp, m.Col, m.t, m.x, rule=MassBalanceSolid_rule)
m.Equilibrium=Constraint(m.Comp, m.Col, m.t, m.x, rule=Equilibrium_rule)

if DV == True:
	def DeadVolume_rule(m, i,j,k):
		if k==0: return Constraint.Skip
		else: return m.dCRedt[i,j,k]/(m.StepTime*m.HT[j]*Nfet) + m.URe[j]*m.dCRedx[i,j,k]/m.LRe == 0
	m.DeadVolume=Constraint(m.Comp, m.t, m.xRe, rule=DeadVolume_rule)

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

if DV == True:
	def CRe_CSS_p_rule(m, i,p):
		return m.CRe[i,1,p] - m.CRe[i,0,p]==0
	m.CRe_CSS_p=Constraint(m.Comp, m.xRe, rule=CRe_CSS_p_rule)

# -------------------------------------------------------------------
# Purity and Recovery Constraints
# -------------------------------------------------------------------

def ExtractRecoveryConstraint_rule(m):
	return sum(m.intCE[2,j] for j in m.Col) >= (m.ExtractRecMin[2]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[2]/m.NCP for j in m.Col for k in m.t)

def ExtractPurityConstraint_rule(m):
	return sum(m.intCE[2,j] for j in m.Col) >= (m.ExtractPurMin[2]/100)*(sum(m.intCE[n,j] for n in m.Comp for j in m.Col))

def RaffinateRecoveryConstraint_rule(m):
	return sum(m.intCR[1,j] for j in m.Col) >= (m.RaffinateRecMin[1]/100)*sum(m.HT[k]*m.UF[j,k]*m.CF[1]/m.NCP for j in m.Col for k in m.t)

def RaffinatePurityConstraint_rule(m):
	return sum(m.intCR[1,j] for j in m.Col) >= (m.RaffinatePurMin[1]/100)*(sum(m.intCR[n,j] for n in m.Comp for j in m.Col))

m.ExtractRecoveryConstraint = Constraint(rule=ExtractRecoveryConstraint_rule)
m.ExtractPurityConstraint= Constraint(rule=ExtractPurityConstraint_rule)
m.RaffinateRecoveryConstraint = Constraint(rule=RaffinateRecoveryConstraint_rule)
m.RaffinatePurityConstraint = Constraint(rule=RaffinatePurityConstraint_rule)