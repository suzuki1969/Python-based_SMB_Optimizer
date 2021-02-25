ObjectiveFunction = '''


@instance.Objective(sense = maximize)
def obj(m):
    return 3600*(sum(m.HT[k]*m.UF[j,k]/m.NCP for j in m.Col for k in m.t if k != 0)) + 1e-5*math.log(value(m.Slack))


Notif = ["----- Throughput Maximizer with Slack Variable is imported -----"]

'''