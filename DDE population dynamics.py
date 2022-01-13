# Delay differential equation model for predicting insect population dynamics under seasonal temperature variation and climate change
# NOTE: If code yields error: "Restarting kernel...", try increasing max_delay in the DDE solver section
# NOTE: If code yields error: "Unsuccessful Integration: Could not integrate with the given tolerance parameters",
#       one of the life history traits is below the minimum tolerance (1e-10)
# NOTE: if code yields error: "CompileError: command 'gcc' failed with exit status 1", one of the parameters is not assigned a value
# NOTE: if code yields error: "IndexError: index 0 is out of bounds for axis 0 with size 0", the location is likely wrong
# NOTE: if code yields warning: "ufunc 'isfinite' not supported for the input types" or "ValueError: cannot convert float NaN to integer",
#       population densities are very high, but this is not a problem for the density-independent model
#       if densities are not plotted, check the output .csv file and remove any rows with 0's after largest density
    
# IMPORT PACKAGES
from numpy import arange, hstack, vstack, savetxt, diff, isnan
from sympy import N
from jitcdde import jitcdde, y, t
from symengine import exp, pi, cos, sin, asin
from matplotlib.pyplot import subplots, xlabel, ylabel, xlim, ylim, yscale, plot, show
from pandas import read_csv
from jitcxde_common import conditional
import os


# SET WORKING DIRECTORY TO SAME AS THE PYTHON CODE
cwd = os.getcwd()
if cwd != '/Users/johnson/Documents/Christopher/GitHub/Johnson_Insect_Responses':
    os.chdir('/Users/johnson/Documents/Christopher/GitHub/Johnson_Insect_Responses')


# USER: Enter species, location, and time period
species = "Hyadaphis pseudobrassicae"
location = "US Columbia"
period = "Historical"
#period = "Future"

# USER: Save data to CSV file?
save_data = False

# USER: Use minimum temperature threshold?
minT = True

# USER: Include competition?
comp = False

# USER: Incorporate diurnal temperature fluctuations?
daily = False

# USER: Is model fit to census data?
census = False

# INPUT TEMPERATURE RESPONSE PARAMETERS AND TEMPERATURE PARAMETERS
data = read_csv("Temperature response parameters.csv")
temp_data = read_csv("Temperature parameters.csv")

# INSECT LIFE HISTORY TRAIT TEMPERATURE RESPONSES
spData = data[data["Species"] == species + " " + location]

# HABITAT TEMPERATURE PARAMETERS
temp_data = temp_data[temp_data["Species"] == species + " " + location]
# diurnal variation (for Tmin)
if period == "Historical":
    diurnal = temp_data["amplD.h"].values[0]
else:
    diurnal = temp_data["amplD.f"].values[0]
if daily == False:
    temp_data = read_csv("Temperature parameters Tave.csv")
    temp_data = temp_data[temp_data["Species"] == species + " " + location]

# DEFINE MODEL PARAMETERS
# Time parameters
yr = 365 # days in year
start_date = 0 # day on which to start model
if comp == True:
    init_years = 0 # how many years into climate change to start model
    max_years = init_years + 75 # how long to run simulations
else:
    init_years = 65 # how many years into climate change to start model
    max_years = 10 # how long to run simulations
tstep = 1 # time step = 1 day
CC_years = max_years # how long before climate change "equilibrates"

# Initial abundances
initJ = 100.
initA = 100.

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

# Life history and competitive traits
# fecundity
bTopt = spData["bTopt"].values[0]
Toptb = spData["Toptb"].values[0]
sb = spData["sb"].values[0]
# maturation
mTR = spData["mTR"].values[0]
TR = spData["TR"].values[0]
AmJ = spData["AmJ"].values[0]
skew = spData["skew"].values[0]
AL = spData["AL"].values[0]
TL = spData["TL"].values[0]
AH = spData["AH"].values[0]
TH = spData["TH"].values[0]
Tmin = spData["Tmin"].values[0] #+ 2*abs(diurnal) # minimum developmental temperature
if census == True: # if model is fit to census data, use temperature at start of experiment as Tmin
    Tmin = temp_data["Tstart"].values[0]
# mortality
dJTR = spData["dJTR"].values[0]
AdJ = spData["AdJ"].values[0]
dATR = spData["dATR"].values[0]
AdA = spData["AdA"].values[0]
# competition
if comp == True:
    qTopt = spData["qTopt"].values[0]
else:
    qTopt = 0
Toptq = Toptb #spData["Toptq"].values[0]
sq = sb #spData["sq"].values[0]
#Aq = spData["Aq"].values[0]
#Tmax = spData["Tmax"].values[0]
#qTopt = qTR*exp(Aq*(1/TR - 1/Tmax))


# FUNCTIONS
# Seasonal temperature variation (K) over time
def T(x):
        return conditional(x, start_date, meanT, # during "pre-history", habitat temperature is constant at its mean
                           conditional(x, CC_years*yr, (meanT + delta_mean*(x+init_years*yr)) - (amplT + delta_ampl*(x+init_years*yr)) * cos(2*pi*(x + shiftT)/yr)  - amplD * cos(2*pi*x), # temperature regime during climate change
                                       (meanT + delta_mean*CC_years*yr) - (amplT + delta_ampl*CC_years*yr) * cos(2*pi*(x + shiftT)/yr)  - amplD * cos(2*pi*x))) # temperature regime after climate change "equilibriates"
'''
# Plot temperature function
xvals = arange(0,1*yr,0.1)
yvals = vstack([(meanT + delta_mean*i) - (amplT + delta_ampl*i) * cos(2*pi*(i + shiftT)/yr)  - amplD * cos(2*pi*i) for i in xvals ])
plot(xvals,yvals)
show()
'''


# Life history functions
# fecundity
def b(x):
    return conditional(bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2)), 1e-5, 1e-5, bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2))) # If b(T) < jitcdde min tolerance, then b(T) = 0

# maturation rates
def mJ(x):
    return conditional(mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + skew * (exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x))))), 1e-5,
                       1e-5, mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + skew * (exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x)))))) # If mJ(T) < jitcdde min tolerance, then mJ(T) = 0

# mortality rates
def dJ(x):
    return dJTR * exp(AdJ * (1/TR - 1/T(x)))
def dA(x):
    return conditional(T(x), Tmin, 1, dATR * exp(AdA * (1/TR - 1/T(x)))) # if temperature < developmental min, Tmin, then dA = 1; otherwise, use dA(T(x))

# density-dependence due to competition
def q(x):
    return qTopt * exp(-(T(x)-Toptq)**2/2/sq**2)


# Minimum developmental temperature
#else:
if minT == True:
    def M(x):
        return conditional(T(x), Tmin, 0, 1) # if temperature < developmental min, Tmin, then development M = 0; otherwise, M = 1
else:
    def M(x):
        return 1


# DDE MODEL
# Define state variables
J,A,S,τ = [y(i) for i in range(4)]

# DDE model
f = {
    J: conditional(t, start_date, 0, M(t)*b(t)*A*exp(-q(t)*A) - M(t-τ)*b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dJ(t)*J), # juveniles: y(0)
    
    A: conditional(t, start_date, 0, M(t-τ)*b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dA(t)*A), # Adults: y(1)
    
    S: conditional(t, start_date, 0, S*(mJ(t)/mJ(t-τ)*dJ(t-τ) - dJ(t))), # Through-stage survivorship: y(2)
    
    τ: conditional(t, start_date, 0, 1 - mJ(t)/mJ(t-τ)), # Developmental time-delay: y(3)
    }


# RUN DDE SOLVER
# Time and initial conditions
times = arange(0, max_years*yr, tstep)
init = [ initJ, initA, exp(-dJ(-1e-3)/mTR), 1./mTR]

# Run DDE solver
DDE = jitcdde(f, max_delay=1e5, verbose=False)
DDE.constant_past(init)
DDE.compile_C(simplify=False, do_cse=False, verbose=True)
DDE.adjust_diff()


# SAVE DATA
# array containing time and state variables
data = vstack([ hstack([time, DDE.integrate(time)]) for time in times ])

# set values below 1e-5 or NAN to 0
data[data < 1e-5] = 0
data[isnan(data)] = 0

# save data to csv 
if save_data == True:
    if comp == True and daily == True:
        filename = 'Time series data/' + period + ' time series ' + spData["Species"].values[0] + '.csv'
        savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
    if comp == True and daily == False:
        filename = 'Time series data Tave/' + period + ' time series ' + spData["Species"].values[0] + '.csv'
        savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
    if comp == False and daily == True:
        filename = 'Time series data DI/' + period + ' time series ' + spData["Species"].values[0] + '.csv'
        savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
    if comp == False and daily == False:
        filename = 'Time series data DI Tave/' + period + ' time series ' + spData["Species"].values[0] + '.csv'
        savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')


# PLOT
fig,ax = subplots()
ax.plot(data[:,0], data[:,1], label='J')
ax.plot(data[:,0], data[:,2], label='A')
#ax.plot(data[:,0], data[:,3], label='S')
#ax.plot(data[:,0], data[:,4], label='τ')
ax.legend(loc='best')
xlabel("time (days)")
ylabel("population density")
yscale("log")
xlim((max_years-max_years)*yr,(max_years-0)*yr)
ylim(0,100000000000000000000000000000)

