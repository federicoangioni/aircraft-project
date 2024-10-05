import math
import numpy as np
import matplotlib.pyplot as plt

MTOW = 43090 #kg
OEW = 24593 #kg
fuel_weight = 8832 #kg
b = 28.08
S = 93.5
c4_sweep = math.radians(17.45)
tr = 0.235
L_fus = 35.50

V_app = 128 # knots
h_app = 0
V_cr = 430 # knots
h_cruise = 30000 # ft

c_r = 5.392

V = 128 # knots
T = 288.15 # K
a = np.sqrt(1.4*287.0*T)
M = (V/1.943845249222) / a
alpha_0l = np.radians(-1.8)
alphas = np.radians(np.linspace(-5, 20))
cl_max = 1.340

dc_cf_20deg = 0.286
dc_cf_50deg = 0.5
cf_c = 0.35
flap_factor = 1.6
y0 = 0.18
y1 = 0.6

C_1 = 0.5

beta = np.sqrt(1-M**2)
eta = 0.95

def QCSweep_to_LESweep(QCSweep, taper_ratio, wing_span, root_chord):
    return math.atan(math.tan(QCSweep) + (1/4) * ((2 * root_chord) / wing_span) * (1 - taper_ratio))

def LESweep_to_HalveCordSweep(LESweep, taper_ratio, wing_span, root_chord):
    return math.atan(math.tan(LESweep) - (1/2) * ((2 * root_chord) / wing_span) * (1 - taper_ratio))

def LESweep_to_hingelineSweep(LESweep, taper_ratio, wing_span, root_chord, c_f_div_c):
    return math.atan(math.tan(LESweep) - (c_f_div_c) * ((2 * root_chord) / wing_span) * (1 - taper_ratio))

def deltaC_l_max_f(dc_cf, cf_c, flap_factor):
    cstar_cf = 1/(cf_c) + dc_cf
    cstar_c = cf_c * cstar_cf
    dC_l_max = cstar_c*flap_factor
    return dC_l_max

def chordaty(yratio,Cr,tr):
    c_y = Cr-Cr*(1-tr)*(yratio)
    return c_y

def area_wf(cy0, cy1, yratio0, yratio1, b):
    b_wf = (yratio1-yratio0)*b
    S_wf = b_wf*((cy0+cy1)/2)
    return S_wf

def deltaC_L_max_f(dC_l_max, S_wf, S, Lambda_hl):
    dC_L_max = 0.9*dC_l_max * (S_wf/S) * math.cos(Lambda_hl)
    return dC_L_max

LE_sweep = np.degrees(QCSweep_to_LESweep(c4_sweep, tr, b, c_r))
HC_sweep = np.degrees(LESweep_to_HalveCordSweep(np.radian(LE_sweep), tr, b, c_r))
Lambda_hl = np.degrees(LESweep_to_hingelineSweep(tr,b,c_r,cf_c))

AR = b**2/S
AR_limit = 4/((C_1+1)*math.cos(np.radians(LE_sweep)))

dclalpha = 2*np.pi*AR/(2+np.sqrt(4+((AR*beta/eta)**2))*(1+((np.tan(HC_sweep))**2)/(beta**2)))

print(f"The dcl/dalpha slope is {dclalpha*np.pi/180}")


# C_l max definition for clean wing
# According to DATCOM zero lift angle is going to be the same


# NACA  64(1)-212
dy = 21.3*0.12

DCLmaxDcl= .82
alpha_delta_CL_Max = 2.4
CL_max = .82*cl_max
alpha_s = CL_max/dclalpha + (alpha_0l) + np.radians(alpha_delta_CL_Max)

cy0 = chordaty(y0,c_r,tr)
cy1 = chordaty(y1,c_r,tr)

Swf = area_wf(cy0, cy1, y0, y1, b)
dC_l_max_f_20deg = deltaC_l_max_f(dc_cf_20deg, cf_c, flap_factor)
dC_L_max_f_20deg = deltaC_L_max_f(dC_l_max_f_20deg, Swf, S, Lambda_hl)

dC_l_max_f_50deg = deltaC_l_max_f(dc_cf_50deg, cf_c, flap_factor)

CL = dclalpha*(alphas-alpha_0l)

print(CL_max, np.degrees(alpha_s))

print(f'The leading edge sweep is: {LE_sweep}')
print(f'The leading edge sweep is: {LE_sweep}')
print(f'The lower limit AR for High AR DATCOM method is: {AR_limit}')
print(f'The leading edge sweep is: {LE_sweep}')
print(f'The CL_max is: {CL_max}')

fig, axs = plt.subplots(figsize=(8, 8))
axs.grid(True)
axs.axhline(linewidth=1, color="k")
axs.axvline(linewidth=1, color="k")
x = np.degrees(alphas)
axs.plot(x[:28+1], CL[:28+1], color='black')
axs.scatter(np.degrees(alpha_s), CL_max, color='red', label= '$C_{L_{max}}$')
print(np.degrees(alpha_s))
#axs.scatter(x[35], CL[35])
#coef = np.array([-0.00838662,  0.23708306, -0.57673129])
#vector = np.array([x**2, x, 1])

axs.plot(x[28:42], - 0.00838662*x[28:42]**2 + 0.23708306*x[28:42] - 0.57673129, color = 'black')
plt.legend()
plt.show()