FeConst = '''


instance.Fe_Const = ConstraintList()

if PowerFeed == 'yes':

    @instance.Constraint(instance.t)
    def no_2portF23(m, k):
        return (m.UF[2,k]*3600)*(m.UF[3,k]*3600) <= 1e-6
    @instance.Constraint(instance.t)
    def no_2portF34(m):
        return (m.UF[3,k]*3600)*(m.UF[4,k]*3600) <= 1e-6
    @instance.Constraint(instance.t)
    def no_2portF42(m):
        return (m.UF[4,k]*3600)*(m.UF[2,k]*3600) <= 1e-6

elif PowerFeed == 'no':

    def no_2portF23(m):
        return (m.UF[2,m.t[2]]*3600)*(m.UF[3,m.t[2]]*3600) <= 1e-6
    def no_2portF34(m):
        return (m.UF[3,m.t[2]]*3600)*(m.UF[4,m.t[2]]*3600) <= 1e-6
    def no_2portF42(m):
        return (m.UF[4,m.t[2]]*3600)*(m.UF[2,m.t[2]]*3600) <= 1e-6

instance.Fe_Const.add(no_2portF23(instance))
instance.Fe_Const.add(no_2portF34(instance))
instance.Fe_Const.add(no_2portF42(instance))


Name = ["Fe_Const"]

Notif = ["----- Feed Constraints ver.0.1 is imported -----"]

'''