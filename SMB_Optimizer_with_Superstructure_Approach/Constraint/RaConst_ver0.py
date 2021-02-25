RaConst = '''


if PowerFeed == 'yes':

    @instance.Constraint(instance.t)
    def no_2portR34(m,k):
        return (m.UR[3,k]*3600)*(m.UR[4,k]*3600) <= 1e-6

elif PowerFeed == 'no':

    @instance.Constraint()
    def no_2portR34(m):
        return (m.UR[3,m.t[2]]*3600)*(m.UR[4,m.t[2]]*3600) <= 1e-6


Name = ["no_2portR34"]

Notif = ["----- Raffinate Constraints ver.0 is imported -----"]

'''