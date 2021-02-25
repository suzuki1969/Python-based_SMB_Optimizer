DeConst = '''


@instance.Constraint()
def UDMaxBound(m):
    m.Slack = Var(initialize=1e-2, within=PositiveReals)
    return 3600*sum(m.HT[k]*m.UD[j,k]/m.NCP for j in m.Col for k in m.t if k != 0) + m.Slack == UB


Name = ["UDMaxBound"]

Notif = ["----- Desorbent Constaint with Slack Variable ver.0 is imported -----"]

'''