import math
import numpy as np


c4_sweep = math.radians(17.45)
lambd = 0.235
b = 28.08
c_r = 5.392
S = 93.5
V = 128 # knots
T = 288.15 # K
a = np.sqrt(1.4*287.0*T)
M = (V/1.943845249222) / a

beta = np.sqrt(1-M**2)
eta = 0.95


def QCSweep_to_LESweep(QCSweep, taper_ratio, wing_span, root_chord):
    return math.atan(math.tan(QCSweep) + (1/4) * ((2 * root_chord) / wing_span) * (1 - taper_ratio))

def QCSweep_to_HalveCordSweep(QCSweep, taper_ratio, wing_span, root_chord):
    return math.atan(math.tan(QCSweep) - (1/3) * ((2 * root_chord) / wing_span) * (1 - taper_ratio))

HC_sweep = QCSweep_to_HalveCordSweep(c4_sweep, lambd, b, c_r)
LE_sweep = np.degrees(QCSweep_to_LESweep(c4_sweep, lambd, b, c_r))
AR = b**2/S

dclalpha = 2*np.pi*AR/(2+np.sqrt(4+((AR*beta/eta)**2))*(1+((np.tan(HC_sweep))**2)/(beta**2)))

print(f"The dcl/dalpha slope is {dclalpha*np.pi/180}")
print(2*np.pi*180/np.pi)