from scipy.optimize import least_squares as olsq
from . import solve


def halfGran(Macid,EMF,Tk,Msamp,AT_cert,XT,KX):
    
    return olsq(lambda Cacid: solve.halfGran(Macid,EMF,Tk,Msamp,Cacid,*XT,*KX,
        suppress_warnings=True)[0] - AT_cert,0.1, method='lm')


def MPH(Macid,EMF,Tk,Msamp,AT_cert,XT,KX):
    
    return olsq(lambda Cacid: \
        solve.MPH(Macid,EMF,Tk,Msamp,Cacid,*XT,KX)['x'][0] - AT_cert,
        0.1, method='lm')


def DAA03(Macid,EMF,Tk,Msamp,AT_cert,XT,KX):
    
    return olsq(lambda Cacid: \
        solve.DAA03(Macid,EMF,Tk,Msamp,Cacid,*XT,KX)['x'][0] - AT_cert,
        0.1, method='lm')
