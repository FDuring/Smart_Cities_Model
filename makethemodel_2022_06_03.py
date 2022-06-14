import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem

#-------  Physical values  ------------------------------------------------
# P-controler gain
Kp = 1e4            # almost perfect controller Kp -> ∞
#Kp = 1e-3           # no controller Kp -> 0

# -------  Geometry ----------------- 
#Old
H = 2.4           # m height of the rooms

#Bath room
l_b = 2           #  m length of the bathoom
w_b = 2           #  m width of the bathoom

Va_b = l_b*w_b*H                # m³ volume of air
ACH = 1                         # air changes per hour 
Va_dot_b = ACH * Va_b / 3600    # m³/s air infiltration
'''Va_dot_b VARIABLE NEVER USED'''

#Main room
l = 5             #  m length of the main room
w = 4             #  m width of the main room 

Va = l*w*H - Va_b           # m³ volume of air in the main room substracted by the bathroom
ACH = 1                     # air changes per hour
Va_dot = ACH * Va / 3600    # m³/s air infiltration

#Window
h_w = 1.4
b_w = 1

#doors
h_d = 2
b_d = 1

# Surface areas [m^2]
Sw = h_w * b_w                      #Surface area window
Sd = h_d * b_d                      #Surface area door
Swallouter =  w * H - Sw            #Surface area outer wall minus window
Swallroombath = l_b * w_b * H - Sd  #Surface room to bath minus door
Swallhallbath = (w - w_b) * H       #Surface hall to bath
Swallroomhall = (w - w_b) * H - Sd  #Surface hall to room minus door

# -------  Thermophyscal properties  ------------------------------------------

air = {'Density': 1.2,                      # kg/m³
       'Specific heat': 1006}               # J/kg.K

""" Incropera et al. (2011) Fundamantals of heat and mass transfer, 7 ed,
    Table A3, concrete (stone mix) p. 993, glass plate p.993, insulation polystyrene extruded (R-12) p.990"""



wall = {'Conductivity': [1.4, 1.4, 1.4, 1.4, 0.04, 1.4, 0.014],      # W/m.K
        'Density': [2500, 2500, 2500, 2500, 30, 2500, 400],            # kg/m³
        'Specific heat': [880, 880, 880, 880, 920, 840, 1900],      # J/kg.K
        'Width': [0.3, 0.1, 0.1, 0.1, 0.15, 0.02, 0.05],        #Glas is set to zero!
        'Surface': [Swallouter, Swallroombath, Swallhallbath, Swallroomhall, Swallouter, Sw, Sd],
        'Slice': [4, 2, 2, 2, 2, 1, 1]}                    # in how many parts is do you wan the material segmented

wall = pd.DataFrame(wall, index=['Wallouter', 'Wallroombath', 'Wallhallbath','Wallroomhall', 'Insulation', 'Window', 'Door'])

#------- Radiative properties -------------------------------------------------

""" concrete EngToolbox Emissivity Coefficient Materials """
ε_wLW = 0.9     # long wave wall emmisivity

""" grey to dark surface EngToolbox, Absorbed Solar Radiation by Surface Color """
α_wSW = 0.2     # absortivity white surface

""" Glass, pyrex EngToolbox Absorbed Solar Radiation by Surface Color """
ε_gLW = 0.9     # long wave glass emmisivity

""" EngToolbox Optical properties of some typical glazing mat Window glass """
τ_gSW = 0.83    # short wave glass transmitance

α_gSW = 0.1     # short wave glass absortivity

σ = 5.67e-8     # W/m².K⁴ Stefan-Bolzmann constant
Fwg = 1 / 5     # view factor wall - glass
Tm = 20 + 273   # mean temp for radiative exchange

# convection coefficients, W/m² K https://de.wikipedia.org/wiki/Wärmeübergangskoeffizient
h = pd.DataFrame([{'in': 4, 'out': 10}]) #indoor and outdoor

#Uwindow = pd.DataFrame([{'in': 2.7, 'out': 2.7}])


#------- Thermal circuit -------------------------------------------------
# Thermal conductances

# Conduction
G_cd = wall['Conductivity'] / wall['Width'] * wall['Surface']

# Convection
G_cv = {'in': [], 'out': []}
for i in wall['Surface']:
    G_cv['in'].append(i * h['in'][0])
    G_cv['out'].append(i * h['out'][0])
G_cv = pd.DataFrame(G_cv, index = ['Wallouter', 'Wallroombath', 'Wallhallbath','Wallroomhall', 'Insulation', 'Window', 'Door'])


# Long-wave radiation exchange
GLW1 = ε_wLW / (1 - ε_wLW) * wall['Surface']['Insulation'] * 4 * σ * Tm**3
GLW2 = Fwg * wall['Surface']['Insulation'] * 4 * σ * Tm**3
GLW3 = ε_gLW / (1 - ε_gLW) * wall['Surface']['Window'] * 4 * σ * Tm**3

# long-wave exg. wall-glass
GLW = 1 / (1 / GLW1 + 1 / GLW2 + 1 / GLW3)

# ventilation & advection
#Gv = Va_dot * air['Density'] * air['Specific heat']

# glass: convection outdoor & conduction & convection indoor
#Ggs = float(1 / (1 / G_cv['out']['Window'] + 1 / (2 * G_cd['Window']) + 1 / (2 * G_cd['Window']) + 1 / G_cv['in']['Window']))

#Window combined value U-Value [W/m²]
G_win = 1

#Door combinded value U-Value [W/m²] 
G_door = 1.5

# Thermal capacities
C = wall['Density'] * wall['Specific heat'] * wall['Surface'] * wall['Width']
C['Air'] = air['Density'] * air['Specific heat'] * Va


#----------------------    our MATRIX / our code     -------------------------
# Incidence matrix
A = np.zeros([16, 10]) #create an matrix with zeros
# adds values for the ones /= 0
A[0, 0] = 1 #position 0,0 and 
A[1, 0], A[1, 1] = -1, 1
A[2, 1], A[2, 2] = -1, 1
A[3, 2], A[3, 3] = -1, 1
A[4, 3], A[4, 4] = -1, 1
A[5, 4], A[5, 5] = -1, 1
A[6, 5] = 1
A[7, 5] = 1
A[8, 5], A[8, 6] = 1,-1
A[9, 6], A[9, 7] = -1, 1
A[10, 5], A[10, 7]  = 1, -1
A[11, 8] = 1
A[12, 6], A[12, 8]  = 1, -1
A[13, 9] = 1
A[14, 5], A[14, 9]  = 1, -1
A[15, 5] = 1

#-- temperatures --
#T0 =   10     #  Outside temp, varying based on weather data
#T_isp = 20     #	Indoor temp. setpoint / radiator	
#T_hall = 15    #	Hall temp.	


#----------------------    End of our MATRIX     -----------------------------
# note: just add values in this  vector!
#Matthias: I think we need to rethink our model because of the G matrix! 
#His makes no sense :( but i think we need more resisors

G = np.diag([
        G_cv['out']['Wallouter'], 
        2 * G_cd['Wallouter'], #If sliced in 4 it should be G_cd / 8 -> Tutorial 2
        2 * G_cd['Wallouter'],
        2 * G_cd['Insulation'],
        2 * G_cd['Insulation'],
        G_cv['in']['Wallouter'],
        G_win, 
        Kp, 
        G_door, #G_cv['in']['Door'] + 2 * G_cd['Door'] + 2 * G_cd['Door'] + G_cv['in']['Door'], # 2 * cond / 2 slices ?
        G_cv['in']['Wallroombath'] + 2 * G_cd['Wallroombath'],
        2 * G_cd['Wallroombath'] + G_cv['in']['Wallroombath'],
        G_cv['in']['Wallhallbath'] + 2 * G_cd['Wallhallbath'], 
        2 * G_cd['Wallhallbath'] + G_cv['in']['Wallhallbath'],   
        G_cv['in']['Wallroomhall'] + 2 * G_cd['Wallroomhall'],   
        2 * G_cd['Wallroomhall'] + G_cv['in']['Wallroomhall'],      
        G_door, #G_cv['in']['Door'] + 2 * G_cd['Door'] + 2 * G_cd['Door'] + G_cv['in']['Door'],
        ])


'''
#G original
G = np.diag(
        [Gw.iloc[0]['out'],
        2 * G_cd['Concrete'],
        2 * G_cd['Concrete'],
        2 * G_cd['Insulation'],
        2 * G_cd['Insulation'],
        GLW,
        Gw.iloc[0]['in'],
        Gg.iloc[0]['in'],
        Ggs,
        2 * G_cd['Glass'],
        Gv,
        Kp])
'''
C = np.diag([0, C['Wallouter'], 0, C['Insulation'], 0, 0, 0,
             C['Wallroombath'], C['Wallhallbath'], C['Wallroomhall']])

#------- More written MATRIXES -----------------------------------------------

#b = np.array([T0, 0, 0, 0, 0, 0, T0, T_isp, 0, 0, 0, T_hall, 0, T_hall, 0, T_hall])
b = np.zeros(16)
b[[0, 6, 7, 11, 14, 15]] = 273.15 + np.array([10, 10, 20, 18, 18, 18])
    
    
# Matthias: We need to define values for our heat sources f
#f = np.array([F0, 0, 0, 0, F1, Q_a, Q_b, 0, 0, 0])
f = np.zeros(10)
f[[0, 4, 5, 6]] = np.array([1000, 500, 500, 500]) # just random values for now

y = np.ones(10)

u = np.hstack([b[np.nonzero(b)], f[np.nonzero(f)]])

# Thermal circuit -> state-space
# ==============================
[As, Bs, Cs, Ds] = dm4bem.tc2ss(A, G, b, C, f, y)

# Test: comparison steady-state of thermal circuit and state-space
ytc = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
yss = (-Cs @ np.linalg.inv(As) @ Bs + Ds) @ u

print(np.array_str(yss, precision=3, suppress_small=True))
print(np.array_str(ytc, precision=3, suppress_small=True))
print(f'Max error in steady-state between thermal circuit and state-space:\
 {max(abs(yss - ytc)):.2e}')


# Dynamic simulation
# ==================
# Thermal circuit -> state-space with 1 for b, f, y
b = np.zeros(16)
b[[0, 6, 7, 11, 14, 15]] = 1

f = np.zeros(10)
f[[0, 4, 5, 6]] = 1

y = np.zeros(10)
y[[5]] = 1

[As, Bs, Cs, Ds] = dm4bem.tc2ss(A, G, b, C, f, y)

# Maximum time-step
dtmax = min(-2. / np.linalg.eig(As)[0])
print(f'Maximum time step: {dtmax:.2f} s')
dt = 5
dt = 360
print(f'Time step: {dt:.2f} s')



# Step response
# -------------
duration = 3600 * 24 * 2       # [s]
# number of steps
n = int(np.floor(duration / dt))

t = np.arange(0, n * dt, dt)    # time

# Vectors of state and input (in time)
n_tC = As.shape[0]              # no of state variables (temps with capacity)
# u = [To To To Tsp Phio Phii Qaux Phia]
# u = [T0, T0, T_sp, T_hall, T_hall, T_hall, f_out, f_wall, f_elect, f_towel]

u = np.zeros([10, n])
u[0:2, :] = np.ones([2, n])
u[3:6, :] = np.ones([3, n])
#fill the matrix with the changing T0, and t Hall

#u = np.hstack([b[np.nonzero(b)], f[np.nonzero(f)]])

# initial values for temperatures obtained by explicit and implicit Euler
temp_exp = np.zeros([n_tC, t.shape[0]])
temp_imp = np.zeros([n_tC, t.shape[0]])

I = np.eye(n_tC)
for k in range(n - 1):
    temp_exp[:, k + 1] = (I + dt * As) @\
        temp_exp[:, k] + dt * Bs @ u[:, k]
    temp_imp[:, k + 1] = np.linalg.inv(I - dt * As) @\
        (temp_imp[:, k] + dt * Bs @ u[:, k])

y_exp = Cs @ temp_exp + Ds @  u
y_imp = Cs @ temp_imp + Ds @  u

fig, axs = plt.subplots(3, 1)
axs[0].plot(t / 3600, y_exp.T, t / 3600, y_imp.T)
axs[0].set(ylabel='$T_i$ [°C]', title='Step input: To = 1°C')


# Simulation with weather data
# ----------------------------
filename = 'FRA_Lyon.074810_IWEC.epw'
start_date = '2000-01-03 12:00:00'
end_date = '2000-01-10 12:00:00'

# Read weather data from Energyplus .epw file
[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data
weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather[(weather.index >= start_date) & (
    weather.index < end_date)]

# Solar radiation on a tilted surface
surface_orientation = {'slope': 90,
                       'azimuth': 0,
                       'latitude': 45}
albedo = 0.2
rad_surf1 = dm4bem.sol_rad_tilt_surf(weather, surface_orientation, albedo) # W/m²
rad_surf1['Φt1'] = rad_surf1.sum(axis=1)

# Interpolate weather data for time step dt
data = pd.concat([weather['temp_air'], rad_surf1['Φt1']], axis=1)
data = data.resample(str(dt) + 'S').interpolate(method='linear')
data = data.rename(columns={'temp_air': 'To'})

# Indoor temperature set-point
data['Ti'] = 20 * np.ones(data.shape[0]) # set point inside

# Indoor auxiliary heat flow rate
data['Qa'] = 0 * np.ones(data.shape[0])

#Outdoor temperature
data['Thall'] = 18 * np.ones(data.shape[0]) 

# time
t = dt * np.arange(data.shape[0])

#Matthias: change u matrix - soral flow depends on orientation and slope of the earth
# u = [T0, T0, T_sp, T_hall, T_hall, T_hall, f_out, f_wall, f_elect, f_towel]
# u = [To To To Tsp Phio Phii Qaux Phia]
u = pd.concat([data['To'], 
               data['To'], 
               data['Ti'] , 
               data['To'], 
               data['To'],
               data['To'],
               α_wSW * wall['Surface']['Wallouter'] * data['Φt1'], #the folowing four lines are f
               τ_gSW * α_wSW * wall['Surface']['Window'] * data['Φt1'],
               data['Qa'], 
               data['Qa']], axis=1)


'''
u = pd.concat([data['To'], data['To'], data['To'], data['Ti'], # vector b
               α_wSW * wall['Surface']['Wallouter'] * data['Φt1'], #the folowing four lines are f
               τ_gSW * α_wSW * wall['Surface']['Window'] * data['Φt1'],
               data['Qa'], 
               α_gSW * wall['Surface']['Window'] * data['Φt1']], axis=1)
'''

# initial values for temperatures
temp_exp = 20 * np.ones([As.shape[0], u.shape[0]])

# integration in time 
I = np.eye(As.shape[0])
for k in range(u.shape[0] - 1):
    temp_exp[:, k + 1] = (I + dt * As) @ temp_exp[:, k]\
        + dt * Bs @ u.iloc[k, :]

# Indoor temperature
y_exp = Cs @ temp_exp + Ds @ u.to_numpy().T
# HVAC heat flow
q_HVAC = Kp * (data['Ti'] - y_exp[0, :])

# plot indoor and outdoor temperature
axs[1].plot(t / 3600, y_exp[0, :], label='$T_{indoor}$')
axs[1].plot(t / 3600, data['To'], label='$T_{outdoor}$')
axs[1].set(xlabel='Time [h]',
           ylabel='Temperatures [°C]',
           title='Simulation for weather')
axs[1].legend(loc='upper right')

# plot total solar radiation and HVAC heat flow
axs[2].plot(t / 3600,  q_HVAC, label='$q_{HVAC}$')
axs[2].plot(t / 3600, data['Φt1'], label='$Φ_{total}$')
axs[2].set(xlabel='Time [h]',
           ylabel='Heat flows [W]')
axs[2].legend(loc='upper right')

fig.tight_layout()
