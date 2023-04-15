from libraries import *


re = 3.87
e = 15.20 * 8.06554     # в см-1
b=6.5
C0=1.02*10**5 * 8.06554
b1=-2.5457689451844066
b2=5.665376678222135
b3=-4.276960450638895
b4=-2.2417153960400404
A2=8.4*10**6 * 8.06554
alpha=3.37
a2 = 0.26
A4 = 1*10**6 * 8.06554
rm = 4.05
x1=1.12
x2=1.55


T = 300 * 1.380649 * 10**(-16) / (1.6 * 10**(-12)) * 8065.54    # в см-1
h=sp.h * 10**7 / (2 * np.pi)
mxe = 131.293 * 1.6605402 * 10**(-24)
mO2 = 32 * 1.6605402 * 10**(-24)
IO2 = mO2 / 2 * (1.20752 * 10**(-8))**2
mu_O2 = mO2 / 2
mu_xeO2 = mO2 * mxe / (mO2 + mxe)

rO2 = 1.2074
vO2 = 1580