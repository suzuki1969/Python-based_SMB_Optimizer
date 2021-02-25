ExConst = '''


if PowerFeed == 'yes':

    @instance.Constraint(instance.t)
    def no_2portE14(m,k):
        return (m.UE[1,k]*3600)*(m.UE[4,k]*3600) <= 1e-6

elif PowerFeed == 'no':

    @instance.Constraint()
    def no_2portE14(m):
        return (m.UE[1,m.t[2]]*3600)*(m.UE[4,m.t[2]]*3600) <= 1e-6


Name = ["no_2portE14"]

Notif = ["----- Extract Constraints ver.0 is imported -----"]

'''