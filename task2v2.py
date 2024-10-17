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
c_t = c_r * tr

V = 128 # knots
T = 288.15 # K
a = np.sqrt(1.4*287.0*T)
M = (V/1.943845249222) / a
alpha_0l = np.radians(-1.8)
alphas = np.radians(np.linspace(-5, 20))
cl_max = 1.340

dc_cf_takeoff = 0.37
dc_cf_landing = 0.64
cf_c = 0.35
flap_factor = 1.6
y0 = 0.18
y1 = 0.6

C_1 = 0.5

beta = 0.981
eta = 0.95

def QCSweep_to_LESweep(QCSweep, taper_ratio, wing_span, root_chord):
    return math.atan(math.tan(QCSweep) + ((1 * root_chord) / (2*wing_span)) * (1 - taper_ratio))

def LESweep_to_HalveCordSweep(LESweep, taper_ratio, wing_span, root_chord):
    return math.atan(math.tan(np.radians(LESweep)) - (root_chord / wing_span) * (1 - taper_ratio))

def LESweep_to_hingelineSweep(LESweep, taper_ratio, wing_span, root_chord, c_f_div_c):
    return math.atan(math.tan(np.radians(LESweep)) - (1-c_f_div_c) * ((2 * root_chord) / wing_span) * (1 - taper_ratio))

def deltaC_l_max_f(dc_cf, cf_c, flap_factor):
    cstar_cf = 1/(cf_c) + dc_cf
    cstar_c = cf_c * cstar_cf
    dC_l_max = cstar_c*flap_factor
    return dC_l_max, cstar_c

def chordaty(yratio,Cr,tr):
    c_y = Cr-Cr*(1-tr)*(yratio)
    return c_y

def area_wf(cy0, cy1, yratio0, yratio1, b):
    b_wf = (yratio1-yratio0)*b
    S_wf = b_wf*((cy0+cy1)/2)
    return S_wf

def deltaC_L_max_f(dC_l_max, S_wf, S, Lambda_hl):
    dC_L_max = 0.9*dC_l_max * (S_wf/S) * math.cos(np.radians(Lambda_hl))
    return dC_L_max

def dalpha_0L(dalpha_0l, S_wf, S, Lambda_hl):
    dC_L_max = dalpha_0l * (S_wf/S) * math.cos(np.radians(Lambda_hl))
    return dC_L_max

def S_sec(c0,c1,b,y0,y1):
    return b*(y1-y0)*(c0+c1)/2

def CL_alpha_flapped(S_flapped, S, C_L_alpha_clean):
    return (S_flapped/S)* C_L_alpha_clean

LE_sweep = np.degrees(QCSweep_to_LESweep(c4_sweep, tr, b, c_r))
HC_sweep = np.degrees(LESweep_to_HalveCordSweep(np.radian(LE_sweep), tr, b, c_r))
Lambda_hl = np.degrees(LESweep_to_hingelineSweep(tr,b,c_r,cf_c))

AR = b**2/S
AR_limit = 4/((C_1+1)*math.cos(np.radians(LE_sweep)))

dclalpha = 2*np.pi*AR/(2+np.sqrt((4+((AR*beta/eta)**2))*(1+((np.tan(np.radians(HC_sweep)))**2)/(beta**2))))

print(f"The dcl/dalpha slope is {np.degrees(dclalpha**-1)**-1}")
print(f"The dcl/dalpha slope is {dclalpha}")

# C_l max definition for clean wing
# According to DATCOM zero lift angle is going to be the same


# NACA  64(1)-212
dy = 21.3*0.12

DCLmaxDcl= .82
alpha_delta_CL_Max = 2.4
CL_max = DCLmaxDcl*cl_max
alpha_s = CL_max/dclalpha + (alpha_0l) + np.radians(alpha_delta_CL_Max)

dalpha_0l_landing = -15 #deg
dalpha_0l_takeoff = -10 #deg

cy0 = chordaty(y0,c_r,tr)
cy1 = chordaty(y1,c_r,tr)

Swf = area_wf(cy0, cy1, y0, y1, b)
dC_l_max_f_takeoff = deltaC_l_max_f(dc_cf_takeoff, cf_c, flap_factor)[0]
cstar_c_takeoff = deltaC_l_max_f(dc_cf_takeoff, cf_c, flap_factor)[1]
#dC_L_max_f_takeoff = deltaC_L_max_f(dC_l_max_f_takeoff, Swf, S, Lambda_hl)

dC_l_max_f_landing = deltaC_l_max_f(dc_cf_landing, cf_c, flap_factor)[0]
cstar_c_landing = deltaC_l_max_f(dc_cf_landing, cf_c, flap_factor)[1]
#dC_L_max_f_landing = deltaC_L_max_f(dC_l_max_f_landing, Swf, S, Lambda_hl)

dalpha_0L_takeoff = dalpha_0L(dalpha_0l_takeoff, Swf, S, Lambda_hl)
dalpha_0L_landing = dalpha_0L(dalpha_0l_landing, Swf, S, Lambda_hl) 

CL = dclalpha*(alphas-alpha_0l)

alpha_0L_land = alpha_0l+np.radians(dalpha_0L_landing)
C_L_max_land = CL_max+dC_L_max_f_landing

alpha_s_land = C_L_max_land/CL_alpha_land + alpha_0L_land + np.radians(alpha_delta_CL_Max)

print(CL_max, np.degrees(alpha_s))

print(f'The leading edge sweep is: {LE_sweep}')
print(f'The Half Chord sweep is: {HC_sweep}')
print(f'The lower limit AR for High AR DATCOM method is: {AR_limit}')
print(f'The CL_max is: {CL_max}')
print(f'The alpha_stall is: {np.degrees(alpha_s)}')
print(f'The Aspect ratio is: {AR}')
print(f'Beta is: {beta}')
"""
print(f'Hingeline sweep is: {Lambda_hl}')
print(f'Flap cr is: {cy0}')
print(f'Flap ct is: {cy1}')
print(f'Flapped wing area is: {Swf}')
print(f'dCL_f_takeoff is: {dC_L_max_f_takeoff}')
print(f'dCL_f_landing is: {dC_L_max_f_landing}')
print(f'Cstar_c_takeoff is: {cstar_c_takeoff}')
print(f'cstar_c_landing is: {cstar_c_landing}')
print(f'dalpha_takeoff is: {dalpha_0L_takeoff}')
print(f'dalpha_landing is: {dalpha_0L_landing}')


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
#plt.show()