FeConst = '''


@m.Constraint()
def no_2portF2(m):
    return m.UF[2,m.t[2]]*(m.UF[3,m.t[2]]+m.UF[4,m.t[2]]) == 0.0
@m.Constraint()
def no_2portF3(m):
    return m.UF[3,m.t[2]]*(m.UF[4,m.t[2]]+m.UF[2,m.t[2]]) == 0.0


Name = ["define_UFsum", "no_2portF2", "no_2portF3"]

Notif = ["----- Feed Constraints ver.2 is imported -----"]

'''