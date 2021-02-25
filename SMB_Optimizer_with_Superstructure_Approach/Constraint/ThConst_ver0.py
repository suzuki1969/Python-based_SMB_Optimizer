ThConst = '''


@instance.Constraint()
def MinThroughput(m):
	return 3600*sum(m.HT[k]*m.UF[j,k]/m.NCP for j in m.Col for k in m.t if k != 0) >= LB


Name = ["MinThroughput"]

Notif = ["----- Throughput Constaint ver.0 is imported -----"]

'''