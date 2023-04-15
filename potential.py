from constants import *


def V01(r):
    x=r/rm
    return e * (np.exp(-2 * b * (x - 1)) - 2 * np.exp(-b * (x - 1)))

def V02(r):
    x=r/rm
    return e * (b1+b2*(x-x1)+(x-x2)*(b3+(x-x1)*b4))

def V03(r):
    return -C0 / r**6

def V2(r):
    return A2 * np.exp(-alpha * r)-a2 * C0 / (r**6)

def V4(r):
    return A4*np.exp(-alpha*r)

def P2(x):
    return 0.5*(3*x**2-1)

def P4(x):
    return 1/8*(35*x**4-30*x**2+3)

def V10(r,th):
    return (V01(r)+V2(r)*P2(th)+V4(r)*P4(th))
def V20(r,th):
    return (V02(r)+V2(r)*P2(th)+V4(r)*P4(th))
def V30(r,th):
    return (V03(r)+V2(r)*P2(th)+V4(r)*P4(th))

def V(r, th):
    if r/rm<=x1:
        return V10(r, th)
    elif r/rm>x2:
        return V30(r, th)
    else:
        return V20(r, th)
