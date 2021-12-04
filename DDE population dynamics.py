# Delay differential equation model for predicting insect population dynamics under seasonal temperature variation and climate change
# NOTE: If code yields error: "Restarting kernel...", try increasing max_delay in the DDE solver section
# NOTE: If code yields error: "Unsuccessful Integration: Could not integrate with the given tolerance parameters",
#       one of the life history traits is below the minimum tolerance (1e-10)
# NOTE: if code yields error: "CompileError: command 'gcc' failed with exit status 1", one of the parameters is not assigned a value
# NOTE: if code yields error: "IndexError: index 0 is out of bounds for axis 0 with size 0", the location is likely wrong
    
# IMPORT PACKAGES
from numpy import arange, hstack, vstack, savetxt
from sympy import N
from jitcdde import jitcdde, y, t
from symengine import exp, pi, cos
from matplotlib.pyplot import subplots, xlabel, ylabel, xlim, ylim, yscale, plot, show
from pandas import read_csv
from jitcxde_common import conditional
import os


# SET WORKING DIRECTORY TO SAME AS THE PYTHON CODE
cwd = os.getcwd()
if cwd != '/Users/johnson/Documents/Christopher/GitHub/Johnson_Insect_Responses':
    os.chdir('/Users/johnson/Documents/Christopher/GitHub/Johnson_Insect_Responses')


# INPUT TEMPERATURE RESPONSE PARAMETERS AND TEMPERATURE PARAMETERS
data = read_csv("Temperature response parameters.csv")
temp_data = read_csv("Temperature parameters.csv")


# ENTER SPECIES, LOCATION, AND TIME PERIOD
species = "Toxoptera citricidus on C. unshiu"
location = "Japan"
period = "Historical"

# USER: Save data to CSV file?
save_data = True

# USER: Model egg stage separately from juvenile stage?
egg = False

# USER: Incorporate resource variation due to precipitation?
res = False

# USER: Use minimum temperature threshold?
minT = True

# USER: Incorporate diurnal temperature fluctuations?
daily = True


# INSECT SPECIES
spData = data[data["Species"] == species + " " + location]

# LOCATION
temp_data = temp_data[temp_data["Species"] == species + " " + location]
#temp_data = temp_data[temp_data["Species"] == species + " Nigeria"]


# DEFINE MODEL PARAMETERS
# Time parameters
yr = 365 # days in year
init_years = 10 # how many years to use for model initiation
max_years = init_years+80 # how long to run simulations
tstep = 1 # time step = 1 day
CC_years = max_years # how long before climate change "equilibrates"

# Initial abundances
initE = 1.
initJ = 1.
initA = 1.
A0 = 0.01 # initial adult density for calculating per capita population growth rate

# Temperature parameters
if period == "Historical":
    meanT = temp_data["meanT.h"].values[0]
    amplT = temp_data["amplT.h"].values[0] 
    shiftT = temp_data["shiftT.h"].values[0]
    delta_mean = temp_data["delta_mean.h"].values[0]
    delta_ampl = temp_data["delta_ampl.h"].values[0]
    amplD = temp_data["amplD.h"].values[0]
else:
    meanT = temp_data["meanT.f"].values[0]
    amplT = temp_data["amplT.f"].values[0] 
    shiftT = temp_data["shiftT.f"].values[0]
    delta_mean = temp_data["delta_mean.f"].values[0]
    delta_ampl = temp_data["delta_ampl.f"].values[0]
    amplD = temp_data["amplD.f"].values[0]
if daily == False:
    amplD = 0

# Resource parameters
if period == "Historical":
    if res == True:
        meanP = temp_data["meanP.h"].values[0]
        amplP = temp_data["amplP.h"].values[0] 
        shiftP = temp_data["shiftP.h"].values[0]
else:
    if res == True:
        meanP = temp_data["meanP.f"].values[0]
        amplP = temp_data["amplP.f"].values[0] 
        shiftP = temp_data["shiftP.f"].values[0]

# Life history and competitive traits
# fecundity
bTopt = spData["bTopt"].values[0]
Toptb = spData["Toptb"].values[0]
sb = spData["sb"].values[0]
# egg maturation
if egg == True:
    mETR = spData["mETR"].values[0]
    AmE = spData["AmE"].values[0]
    ALE = spData["ALE"].values[0]
    TLE = spData["TLE"].values[0]
    AHE = spData["AHE"].values[0]
    THE = spData["THE"].values[0]
# maturation
mTR = spData["mTR"].values[0]
TR = spData["TR"].values[0]
AmJ = spData["AmJ"].values[0]
skew = spData["skew"].values[0]
AL = spData["AL"].values[0]
TL = spData["TL"].values[0]
AH = spData["AH"].values[0]
TH = spData["TH"].values[0]
Tmin = spData["Tmin"].values[0] # minimum developmental temperature
# mortality
if egg == True:
    dETR = spData["dETR"].values[0]
    AdE = spData["AdE"].values[0]
dJTR = spData["dJTR"].values[0]
AdJ = spData["AdJ"].values[0]
dATR = spData["dATR"].values[0]
AdA = spData["AdA"].values[0]
# competition
qTopt = 0.1*spData["qTopt"].values[0]
Toptq = Toptb #spData["Toptq"].values[0]
sq = sb #spData["sq"].values[0]
#Aq = spData["Aq"].values[0]
#Tmax = spData["Tmax"].values[0]
#qTopt = qTR*exp(Aq*(1/TR - 1/Tmax))


# FUNCTIONS
# Seasonal temperature variation (K) over time
def T(x):
        return conditional(x, 0, meanT, # during "pre-history" (t<0), habitat temperature is constant at its mean
                       conditional(x, init_years*yr, meanT - amplT * cos(2*pi*(x + shiftT)/yr) - amplD * cos(2*pi*x), # during model initiation, delta_mean = 0 and delta_ampl = 0
                                   conditional(x, CC_years*yr, (meanT + delta_mean*(x-init_years*yr)) - (amplT + delta_ampl*(x-init_years*yr)) * cos(2*pi*((x-init_years*yr) + shiftT)/yr)  - amplD * cos(2*pi*(x-init_years*yr)), # temperature regime during climate change
                                               (meanT + delta_mean*CC_years*yr) - (amplT + delta_ampl*CC_years*yr) * cos(2*pi*((x-init_years*yr) + shiftT)/yr)  - amplD * cos(2*pi*(x-init_years*yr))))) # temperature regime after climate change "equilibriates"
'''
# Plot temperature function
xvals = arange(0,1*yr,0.1)
yvals = vstack([(meanT + delta_mean*i) - (amplT + delta_ampl*i) * cos(2*pi*(i + shiftT)/yr)  - amplD * cos(2*pi*i) for i in xvals ])
plot(xvals,yvals)
show()
'''

# Seasonal resource variation (R) due to precipitation
def R(x):
        return conditional(x, init_years*yr, 1, # no resource variation during model initiation
                                   conditional(x, CC_years*yr, # temperature regime during climate change
                                               conditional(meanP - amplP * cos(2*pi*((x-init_years*yr) + shiftP)/yr), 0, 0, 1), 1)) # wet season if monthly precipitation > 0 (orinigally used 6 as defined by Köppen climate classification system)
'''
# Plot resource function
xvals = arange(0,1*yr,1)
#yvals = vstack([meanP - amplP * cos(2*pi*(i + shiftP)/yr) for i in xvals ])
yvals = vstack([conditional(meanP - amplP * cos(2*pi*(i + shiftP)/yr), 0, 0, 1) for i in xvals ])
plot(xvals,yvals)
show()
'''

# Life history functions
# fecundity
def b(x):
    if res==True:
        return conditional(R(x) * bTopt * exp(-(T(x)-Toptb)**2/2/sb**2), 1e-5, 1e-5, R(x) * bTopt * exp(-(T(x)-Toptb)**2/2/sb**2)) # If b(T) < jitcdde min tolerance, then b(T) = 0
    else:
        return conditional(bTopt * exp(-(T(x)-Toptb)**2/2/sb**2), 1e-5, 1e-5, bTopt * exp(-(T(x)-Toptb)**2/2/sb**2)) # If b(T) < jitcdde min tolerance, then b(T) = 0

# egg maturation rates
if egg == True:
    def mE(x):
        return conditional(mETR * T(x)/TR * exp(AmE * (1/TR - 1/T(x))) / (1 + skew * (exp(ALE*(1/TLE-1/T(x)))+exp(AHE*(1/THE-1/T(x))))), 1e-5,
                           1e-5, mETR * T(x)/TR * exp(AmE * (1/TR - 1/T(x))) / (1 + skew * (exp(ALE*(1/TLE-1/T(x)))+exp(AHE*(1/THE-1/T(x)))))) # If mJ(T) < jitcdde min tolerance, then mJ(T) = 0

# maturation rates
def mJ(x):
    return conditional(mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + skew * (exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x))))), 1e-5,
                       1e-5, mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + skew * (exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x)))))) # If mJ(T) < jitcdde min tolerance, then mJ(T) = 0

# mortality rates
if egg == True:
    def dE(x):
        return dETR * exp(AdE * (1/TR - 1/T(x)))
def dJ(x):
    return dJTR * exp(AdJ * (1/TR - 1/T(x)))
def dA(x):
    return conditional(T(x), Tmin, 0.5, dATR * exp(AdA * (1/TR - 1/T(x)))) # if temperature < developmental min, Tmin, then dA = dATR; otherwise, use dA(T(x))

# density-dependence due to competition
def q(x):
    return qTopt * exp(-(T(x)-Toptq)**2/2/sq**2)

# Other functions    
# Allee effect
A_thr = 0.00*qTopt # Allee threshold
def Allee(x):
    return conditional(x, A_thr, 0, 1) # if A < A_thr then Allee = 0; otherwise, Allee = 1 

# Minimum developmental temperature
if minT == True:
    def M(x):
        return conditional(T(x), Tmin, 0, 1) # if temperature < developmental min, Tmin, then development M = 0; otherwise, M = 1
else:
    def M(x):
        return 1

# DDE MODEL
# Define state variables
J,A,S,τ,r = [y(i) for i in range(5)]

# DDE model
f = {
    J:  M(t)*b(t)*A*Allee(A)*exp(-q(t)*A) - M(t)*M(t-τ)*b(t-τ)*y(1,t-τ)*Allee(y(1,t-τ))*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dJ(t)*J, # juveniles: y(0)
    
    A:  M(t)*M(t-τ)*b(t-τ)*y(1,t-τ)*Allee(y(1,t-τ))*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dA(t)*A, # Adults: y(1)
    
    S: S*(mJ(t)/mJ(t-τ)*dJ(t-τ) - dJ(t)), # Through-stage survivorship: y(2)
    
    τ: 1 - mJ(t)/mJ(t-τ), # Developmental time-delay: y(3)
    
    r: M(t)*b(t-τ)*A0*Allee(A0)*exp(-q(t-τ)*A0)*mJ(t)/mJ(t-τ)*S - dA(t)*A0 # Low density population growth rate: y(4)
    }


# MODEL WITH EGG STAGE
if egg == True:
    # Define state variables
    E,J,A,SE,SJ,SA,τE,τJ,τA,r = [y(i) for i in range(10)]
    
    # Define functions used in DDE model
    # adult fecundity
    def fA(t):
        return M(t)*b(t)*A*Allee(A)*exp(-q(t)*A)
    # egg development                   
    def gE(t):
        return M(t)*M(t-τE)*b(t-τE)*y(2,t-τE)*Allee(y(2,t-τE))*exp(-q(t-τE)*y(2,t-τE))*mE(t)/mE(t-τE)*SE
    # juvenile development
    def gJ(t):
        return M(t)*M(t-y(6,t-τJ)-τJ)*b(t-y(6,t-τJ)-τJ)*y(2,t-y(6,t-τJ)-τJ)*Allee(y(2,t-y(6,t-τJ)-τJ))*exp(-q(t-y(6,t-τJ)-τJ)*y(2,t-y(6,t-τJ)-τJ))*mE(t-τJ)/mE(t-y(6,t-τJ)-τJ)*mJ(t)/mJ(t-τJ)*y(3,t-τJ)*SJ
    # adult senescence
    def gA(t):
        return 0*M(t)*M(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * b(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * y(2,t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * Allee(y(2,t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA)) * exp(-q(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * y(2,t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA)) * mE(t-y(7,t-τA)-τA)/mE(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * mJ(t-τA)/mJ(t-y(7,t-τA)-τA) * dA(t)/dA(t-τA) * y(3,t-y(7,t-τA)-τA) * y(4,t-τA) * SA

    # DDE model
    f = {
        E: fA(t) - gE(t) - dE(t)*E, # eggs: y(0)

        J: gE(t) - gJ(t) - dJ(t)*J, # juveniles: y(1)
    
        A: gJ(t) - gA(t) - dA(t)*A, # Adults: y(2)
    
        SE: SE*(mE(t)/mE(t-τE)*dE(t-τE) - dE(t)), # Through-egg stage survivorship: y(3)
    
        SJ: SJ*(mJ(t)/mJ(t-τJ)*dJ(t-τJ) - dJ(t)), # Through-juvenile stage survivorship: y(4)
        
        SA: SA*(dA(t)/dA(t-τA)*dA(t-τA) - dA(t)), # Through-adult stage survivorship: y(5)
        
        τE: 1 - mE(t)/mE(t-τE), # Egg developmental time-delay: y(6)
        
        τJ: 1 - mJ(t)/mJ(t-τJ), # Juvenile developmental time-delay: y(7)
        
        τA: 1 - dA(t)/dA(t-τA), # Adult longevity time-delay: y(8)
        
        # Low density population growth rate: y(9)
        r: M(t)*M(t-y(6,t-τJ)-τJ)*b(t-y(6,t-τJ)-τJ)*A0*Allee(A0)*exp(-q(t-y(6,t-τJ)-τJ)*A0)*mE(t-τJ)/mE(t-y(6,t-τJ)-τJ)*mJ(t)/mJ(t-τJ)*y(3,t-τJ)*SJ - dA(t)*A0  -  M(t)*M(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * b(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * A0 * Allee(A0) * exp(-q(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * A0) * mE(t-y(7,t-τA)-τA)/mE(t-y(6,t-y(7,t-τA)-τA)-y(7,t-τA)-τA) * mJ(t-τA)/mJ(t-y(7,t-τA)-τA) * dA(t)/dA(t-τA) * y(3,t-y(7,t-τA)-τA) * y(4,t-τA) * SA
        }


# RUN DDE SOLVER
# Time and initial conditions
times = arange(0, max_years*yr, tstep)
init = [ initJ, initA, exp(-dJ(-1e-3)/mTR), 1./mTR, 0. ]
if egg == True:
    init = [ initE, initJ, initA, exp(-dE(0)/mE(0)), exp(-dJ(0)/mJ(0)), exp(-1), 1./mE(0), 1./mJ(0), 1./dA(0), 0. ]

# Run DDE solver
DDE = jitcdde(f, max_delay=1e5, verbose=False)
DDE.constant_past(init)
DDE.compile_C(simplify=False, do_cse=False, verbose=True)
DDE.adjust_diff()


# array containing time and state variables
data = vstack([ hstack([time, DDE.integrate(time)]) for time in times ])
if egg == False:
    data[:,5] = data[:,5]/data[:,0] # r column from DDE.integrate is actually r*t, so divide by t
    data[0,5] = 0 # reset initial r to 0
else:
    data[:,8] = data[:,8]/data[:,0] # r column from DDE.integrate is actually r*t, so divide by t
    data[0,8] = 0 # reset initial r to 0
    
# SAVE DATA
if save_data == True:
    if egg == False:
        filename = 'Time series data/' + period + ' time series ' + spData["Species"].values[0] + '.csv'
        savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau,r", comments='')
    else:
        #filename = 'Time series data/' + period + ' time series ' + spData["Species"].values[0] + '.csv'
        filename = 'Time series data/' + period + ' time series ' + species + ' Nigeria (egg).csv'
        savetxt(filename, data, fmt='%s', delimiter=",", header="Time,E,J,A,SE,SJ,tauE,tauJ,r", comments='')



# PLOT
fig,ax = subplots()
if egg == False:
    ax.plot(data[:,0], data[:,1], label='J')
    ax.plot(data[:,0], data[:,2], label='A')
    #ax.plot(data[:,0], data[:,3], label='S')
    #ax.plot(data[:,0], data[:,4], label='τ')
else:
    ax.plot(data[:,0], data[:,1], label='E')
    ax.plot(data[:,0], data[:,2], label='J')
    ax.plot(data[:,0], data[:,3], label='A')
ax.legend(loc='best')
xlabel("time (days)")
ylabel("population density")
yscale("linear")
xlim((max_years-1)*yr,(max_years-0)*yr)
ylim(0,1000)
