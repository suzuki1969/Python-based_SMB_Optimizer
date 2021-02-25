ConcentrationConst = '''


@instance.Constraint()
def RaConcentrationConst(m):
	return sum(m.intCR[1,j] for j in m.Col ) >= (m.ConcentrationMin[1]/100)*m.CF[1]*(sum(m.HT[k]*m.UR[j,k]/m.NCP for j in m.Col for k in m.t if k != 0))
@instance.Constraint()
def ExConcentrationConst(m):
	return sum(m.intCE[2,j] for j in m.Col ) >= (m.ConcentrationMin[2]/100)*m.CF[2]*(sum(m.HT[k]*m.UE[j,k]/m.NCP for j in m.Col for k in m.t if k != 0))


Name = ["RaConcentrationConst", "ExConcentrationConst"]

Notif = ["----- Concentration Constraint ver.0 is imported -----"]

'''