# Originally created by Suriya Arulselvan in Dec 2015
# Developed by Zeyu Yang in Jul 2016
# Addtional developement is given by Kensuke Suzuki in 2020
# The mathematical model performed in this file is the same as study by Kawajiri and Biegler
# Bibliographic information: Kawajiri Y, Biegler LT, Optimization Strategies for Simulated Moving Bed and PowerFeed Processes. AIChE Journal. 2006;52:1343-1350.
#                          : Zavala V, Biegler LT, Large-scale parameter estimation in low-density polyethylene tubular reactors. Industrial and Engineering Chemistry Research. 2006;45,23;7867-7881
#                          : Tie et al., Experimental evaluation of simulated moving bed reactor for transesterification reaction synthesis of glycol ether ester. Adsorption. 2019;25,4;795-807

# ==================================================================
# Generalized SMB optimization
# ==================================================================

from pyomo.environ import *
from pyomo.dae import *
from Initdata_SMB_ParamEst import *
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
m.Col = Set(initialize = Col_init, ordered = True)
m.FEX = Set(initialize = FEX_init, ordered = True)
m.CP = Set(initialize = CP_init, ordered = True)
m.Comp = Set(initialize = Comp_init, ordered = True)
m.Data = Set(initialize = Data_init, ordered = True)
m.Regularization = Set(initialize = [1,2])  

# -------------------------------------------------------------------
# Physical Paramters
# -------------------------------------------------------------------

m.L = Param(initialize = ColL)

m.eb = Param(initialize = eb) # void fraction
m.CF = Param(m.Data, m.Comp, initialize = CFeed) # feed concentration
m.H = Var(m.Comp, within=PositiveReals) # Henry's constant
m.Kap = Var(m.Comp, within=PositiveReals) # overall mass transfer coefficient

if Isotype == "Henry":
	pass
elif Isotype == "Langmuir" or Isotype == "anti-Langmuir":
	m.b = Var(m.Comp, within=PositiveReals) # affinity parameter
else:
	sys.exit("ERROR: Confirm Isotype definition")

if Axial_D == True:
	m.Dax = Var(m.Comp, within=PositiveReals) # Axial Dispersion Coefficient
elif Axial_D == False:
	pass
else:
	sys.exit("ERROR: Confirm Axial_D definition")

m.averageCE = Var(m.Data, m.Comp) # average conc in Extract (model)
m.averageCR = Var(m.Data, m.Comp) # average conc in Raffinate (model)

m.CbarE = Param(m.Data, m.Comp, initialize = CExtract) # average conc in Extract (experiment)
m.CbarR = Param(m.Data, m.Comp, initialize = CRaffinate) # average conc in Raffinate (experiment)

m.rho = Param(m.Regularization, initialize = Rho) # regularization parameters  

m.e = Var(m.Data, within=NonNegativeReals)

# -------------------------------------------------------------------
# State Variables
# -------------------------------------------------------------------

m.x = ContinuousSet(bounds=(0,m.L)) # space

m.t = ContinuousSet(bounds=(0,1)) # time

m.HT = Var(m.Data, m.t, bounds = ((0.5/m.NFET), (2/m.NFET))) # length of a finite element in time

m.C = Var(m.Data, m.Comp, m.Col, m.t, m.x) # concentration in liquid phase
m.Q = Var(m.Data, m.Comp, m.Col, m.t, m.x) # concentration in solid phase

if Phase == 'Liquid':
	m.CEQ = Var(m.Data, m.Comp, m.Col, m.t, m.x) # equilibriun concentration in liquid phase
elif Phase == 'Solid':
	m.QEQ = Var(m.Data, m.Comp, m.Col, m.t, m.x) # equilibriun concentration in sloid phase
else:
	sys.exit("ERROR: Confirm Phase definition")

m.intCR = Var(m.Data, m.Comp, m.Col) # integrated concentration of Raffinate
m.intCE = Var(m.Data, m.Comp, m.Col) # integrated concentration of Extract

if DV == True:
	m.xRe = ContinuousSet(bounds=(0,1)) # m.x in dead volume
	m.CRe = Var(m.Data, m.Comp, m.t, m.xRe) # m.C in dead volume
elif DV == False:
    pass
else:
	sys.exit("ERROR: Confirm DeadVolume definition")

# -------------------------------------------------------------------
# Control Variables
# -------------------------------------------------------------------

m.U = Var(m.Data, m.Col, m.t, within=PositiveReals, initialize=1.0) # fluid velocity in column
m.UF = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Feed
m.UD = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Desorbent
m.UR = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Raffinate
m.UE = Var(m.Data, m.Col, m.t, within=PositiveReals) # fluid velocity of Extract

if DV == True:
	m.URe = Var(m.Data, m.t, within=PositiveReals) # fluid velocity in dead volume
	m.LRe = Var(initialize = 0.05*Colum*ColL, within=PositiveReals) # length of dead volume

m.StepTime = Var(m.Data, bounds = (5,3600)) # step time

# -------------------------------------------------------------------
# State Derivatives
# -------------------------------------------------------------------

m.dCdt = DerivativeVar(m.C, wrt=m.t) # differentiated liquid concentration respect to time
m.dQdt = DerivativeVar(m.Q, wrt=m.t) # differentiated solid concentration respect to time
m.dCdx = DerivativeVar(m.C, wrt=m.x) # differentiated liquid concentration respect to space

if Axial_D == True:
	m.dCdx2 = DerivativeVar(m.C, wrt=(m.x, m.x)) # second-order differentiated liquid concentration with respect to space

if DV == True:
	m.dCRedt = DerivativeVar(m.CRe, wrt=m.t)
	m.dCRedx = DerivativeVar(m.CRe, wrt=m.xRe)

# -------------------------------------------------------------------
# Flow rate balances
# -------------------------------------------------------------------

def FlowBalance_rule(m,d,j,k):
	if j==m.NColumn: return Constraint.Skip
	else: return m.U[d,j+1,k]==m.U[d,j,k]-m.UR[d,j,k]-m.UE[d,j,k]+m.UD[d,j+1,k]+m.UF[d,j+1,k]

m.FlowBalance=Constraint(m.Data, m.Col, m.t, rule=FlowBalance_rule)

if DV == True:
	def FlowBalanceD1_rule(m,d,k):
		return m.U[d,1,k] == m.URe[d,k] + m.UD[d,1,k] + m.UF[d,1,k]
	m.FlowBalanceD1=Constraint(m.Data, m.t, rule=FlowBalanceD1_rule)

	def FlowBalance4D_rule(m,d,k):
		return m.URe[d,k] == m.U[d,Colum,k] - m.UR[d,Colum,k] - m.UE[d,Colum,k]
	m.FlowBalance4D=Constraint(m.Data, m.t, rule=FlowBalance4D_rule)
elif DV == False:
	def FlowBalance4D_rule(m,d,k):
		return m.U[d,1,k]==m.U[d,Colum,k]-m.UR[d,Colum,k]-m.UE[d,Colum,k]+m.UD[d,1,k]+m.UF[d,1,k]
	m.FlowBalance4D=Constraint(m.Data, m.t, rule=FlowBalance4D_rule)

# -------------------------------------------------------------------
# Boundary Conditions
# -------------------------------------------------------------------

def BoundaryCondition_rule(m, d,i,j,k):
   if k == 0.0: return Constraint.Skip
   if j==m.NColumn: return Constraint.Skip
   else: return m.C[d,i,j+1,k,0]*m.U[d,j+1,k]==m.C[d,i,j,k,m.L]*(m.U[d,j,k] - m.UE[d,j,k] - m.UR[d,j,k]) + m.CF[d,i]*m.UF[d,j+1,k]

m.BoundaryCondition=Constraint(m.Data, m.Comp, m.Col, m.t, rule=BoundaryCondition_rule)


if DV == True:
	def BoundaryCondition4D_rule(m,d,i,k):
		if k == 0 : return Constraint.Skip 
		return m.CRe[d,i,k,0]*m.URe[d,k]==m.C[d,i,Colum,k,m.L]*(m.U[d,Colum,k] - m.UE[d,Colum,k] - m.UR[d,Colum,k])
	m.BoundaryCondition4D=Constraint(m.Data, m.Comp, m.t, rule=BoundaryCondition4D_rule)

	def BoundaryConditionD1_rule(m,d,i,k):
		if k == 0 : return Constraint.Skip 
		return m.C[d,i,1,k,0]*m.U[d,1,k]==m.CRe[d,i,k,1]*m.URe[d,k] + m.CF[d,i]*m.UF[d,1,k]
	m.BoundaryConditionD1=Constraint(m.Data, m.Comp, m.t, rule=BoundaryConditionD1_rule)
elif DV == False:
	def BoundaryCondition1_rule(m,d,i,k):
		if k == 0.0: return Constraint.Skip
		return m.C[d,i,1,k,0]*m.U[d,1,k]==m.C[d,i,Colum,k,m.L]*(m.U[d,Colum,k]-m.UE[d,Colum,k]-m.UR[d,Colum,k]) + m.CF[d,i]*m.UF[d,1,k]
	m.BoundaryCondition1=Constraint(m.Data, m.Comp, m.t, rule=BoundaryCondition1_rule)

# -------------------------------------------------------------------
# Model equations
# -------------------------------------------------------------------

if Axial_D == True:
	def MassBalanceLiquid_rule(m, d,i,j,k,l):
		if l==0: return Constraint.Skip
		else: return m.eb*m.dCdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*m.NFET) + (1-m.eb)*m.dQdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*m.NFET) + m.U[d,j,k]*m.dCdx[d,i,j,k,l] == m.Dax[i]*m.dCdx2[d,i,j,k,l] 
elif Axial_D == False:
	def MassBalanceLiquid_rule(m, d,i,j,k,l):
		if l==0: return Constraint.Skip
		else: return m.eb*m.dCdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*m.NFET) + (1-m.eb)*m.dQdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*m.NFET) + m.U[d,j,k]*m.dCdx[d,i,j,k,l] == 0 

if Phase == 'Liquid':
	def MassBalanceSolid_rule(m, d,i,j,k,l):
		return (1.0-m.eb)*m.dQdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*m.NFET) == m.Kap[i]*(m.C[d,i,j,k,l]-m.CEQ[d,i,j,k,l])

	if Isotype == "Henry":
		def Equilibrium_rule(m, d,i,j,k,l):
			return m.Q[d,i,j,k,l] == m.H[i]*m.CEQ[d,i,j,k,l]
	elif Isotype == "Langmuir":
		def Equilibrium_rule(m, d,i,j,k,l):
			return m.Q[d,i,j,k,l]*(1.0 + m.b[1]*m.CEQ[d,1,j,k,l] + m.b[2]*m.CEQ[d,2,j,k,l]) == m.H[i]*m.CEQ[d,i,j,k,l]
	elif Isotype == "anti-Langmuir":
		def Equilibrium_rule(m, d,i,j,k,l):
			return m.Q[d,i,j,k,l]*(1.0 - m.b[1]*m.CEQ[d,1,j,k,l] - m.b[2]*m.CEQ[d,2,j,k,l]) == m.H[i]*m.CEQ[d,i,j,k,l]

elif Phase == 'Solid':
	def MassBalanceSolid_rule(m, d,i,j,k,l):
		return m.dQdt[d,i,j,k,l]/(m.StepTime[d]*m.HT[d,k]*m.NFET) == m.Kap[i]*(m.QEQ[d,i,j,k,l]-m.Q[d,i,j,k,l])

	if Isotype == "Henry":
		def Equilibrium_rule(m, d,i,j,k,l):
			return m.QEQ[d,i,j,k,l] == m.H[i]*m.C[d,i,j,k,l]
	elif Isotype == "Langmuir":
		def Equilibrium_rule(m, d,i,j,k,l):
			return m.QEQ[d,i,j,k,l]*(1.0 + m.b[1]*m.C[d,1,j,k,l] + m.b[2]*m.C[d,2,j,k,l]) == m.H[i]*m.C[d,i,j,k,l]
	elif Isotype == "anti-Langmuir":
		def Equilibrium_rule(m, d,i,j,k,l):
			return m.QEQ[d,i,j,k,l]*(1.0 - m.b[1]*m.C[d,1,j,k,l] - m.b[2]*m.C[d,2,j,k,l]) == m.H[i]*m.C[d,i,j,k,l]

m.MassBalanceLiquid=Constraint(m.Data, m.Comp, m.Col, m.t, m.x, rule=MassBalanceLiquid_rule)
m.MassBalanceSolid=Constraint(m.Data, m.Comp, m.Col, m.t, m.x, rule=MassBalanceSolid_rule)
m.Equilibrium=Constraint(m.Data, m.Comp, m.Col, m.t, m.x, rule=Equilibrium_rule)

if DV == True:
	def DeadVolume_rule(m, d,i,j,k):
		if k==0: return Constraint.Skip
		else: return m.dCRedt[d,i,j,k]/(m.StepTime[d]*m.HT[d,j]*m.NFET) + m.URe[d,j]*m.dCRedx[d,i,j,k]/m.LRe == 0
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
	else: return m.C[d,i,j,1,p] - m.C[d,i,j-1,0,p] == m.e[d]

def Q_CSS_p_rule(m, d,i,j,p):
	if j==1: return Constraint.Skip
	else: return m.Q[d,i,j,1,p] - m.Q[d,i,j-1,0,p] == m.e[d]

def C_CSS_1_p_rule(m, d,i,p):
	return m.C[d,i,1,1,p] - m.C[d,i,Colum,0,p]== m.e[d]

def Q_CSS_1_p_rule(m, d,i,p):
	return m.Q[d,i,1,1,p] - m.Q[d,i,Colum,0,p]== m.e[d]

m.C_CSS_p=Constraint(m.Data, m.Comp, m.Col, m.x, rule=C_CSS_p_rule)
m.Q_CSS_p=Constraint(m.Data, m.Comp, m.Col, m.x, rule=Q_CSS_p_rule)
m.C_CSS_1_p=Constraint(m.Data, m.Comp, m.x, rule=C_CSS_1_p_rule)
m.Q_CSS_1_p=Constraint(m.Data, m.Comp, m.x, rule=Q_CSS_1_p_rule)

if DV == True:
	def CRe_CSS_p_rule(m, d,i,p):
		return m.CRe[d,i,1,p] - m.CRe[d,i,0,p]== m.e[d]
	m.CRe_CSS_p=Constraint(m.Data, m.Comp, m.xRe, rule=CRe_CSS_p_rule)
