FeConst = '''


@m.Constraint()
def define_UFsum(m):
    m.UFsum = Var(within=PositiveReals)
    return m.UFsum == m.UF[3,m.t[2]] + m.UF[4,m.t[2]]
@m.Constraint()
def define_UFsum(m):
    m.UFsum = Var(within=PositiveReals)
    return m.UFsum == m.UF[3,m.t[2]] + m.UF[4,m.t[2]]
@m.Constraint()
def no_2portF4(m):
    return m.UFsum*m.UF[4,m.t[2]] == m.UFsum**2


Name = ["define_UFsum", "no_2portF3", "no_2portF4"]

Notif = ["----- Feed Constraints ver.1 is imported -----"]

'''