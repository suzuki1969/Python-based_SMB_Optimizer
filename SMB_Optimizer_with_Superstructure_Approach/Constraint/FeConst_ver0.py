FeConst = '''


if PowerFeed == 'yes':

    @instance.Constraint(instance.t)
    def no_2portF12(m,k):
        return (m.UF[1,k]*3600)*(m.UF[2,k]*3600) <= 1e-6
    @instance.Constraint(instance.t)
    def no_2portF13(m,k):
        return (m.UF[1,k]*3600)*(m.UF[3,k]*3600) <= 1e-6
    @instance.Constraint(instance.t)
    def no_2portF14(m,k):
        return (m.UF[1,k]*3600)*(m.UF[4,k]*3600) <= 1e-6
    @instance.Constraint(instance.t)
    def no_2portF23(m,k):
        return (m.UF[2,k]*3600)*(m.UF[3,k]*3600) <= 1e-6
    @instance.Constraint(instance.t)
    def no_2portF24(m,k):
        return (m.UF[2,k]*3600)*(m.UF[4,k]*3600) <= 1e-6
    @instance.Constraint(instance.t)
    def no_2portF34(m,k):
        return (m.UF[3,k]*3600)*(m.UF[4,k]*3600) <= 1e-6

elif PowerFeed == 'no':

    @instance.Constraint()
    def no_2portF12(m):
        return (m.UF[1,m.t[2]]*3600)*(m.UF[2,m.t[2]]*3600) <= 1e-6
    @instance.Constraint()
    def no_2portF13(m):
        return (m.UF[1,m.t[2]]*3600)*(m.UF[3,m.t[2]]*3600) <= 1e-6
    @instance.Constraint()
    def no_2portF14(m):
        return (m.UF[1,m.t[2]]*3600)*(m.UF[4,m.t[2]]*3600) <= 1e-6
    @instance.Constraint()
    def no_2portF23(m):
        return (m.UF[2,m.t[2]]*3600)*(m.UF[3,m.t[2]]*3600) <= 1e-6
    @instance.Constraint()
    def no_2portF24(m):
        return (m.UF[2,m.t[2]]*3600)*(m.UF[4,m.t[2]]*3600) <= 1e-6
    @instance.Constraint()
    def no_2portF34(m):
        return (m.UF[3,m.t[2]]*3600)*(m.UF[4,m.t[2]]*3600) <= 1e-6


Name = ["no_2portF12", "no_2portF13", "no_2portF14", "no_2portF23", "no_2portF24", "no_2portF34"]

Notif = ["----- Feed Constraints ver.0 is imported -----"]

'''