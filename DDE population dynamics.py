# Delay differential equation model for predicting insect population dynamics under seasonal temperature variation and climate change
# NOTE: If code yields error: "Restarting kernel...", try increasing max_delay in the DDE solver section
# NOTE: If code yields error: "Unsuccessful Integration: Could not integrate with the given tolerance parameters",
#       one of the life history traits is below the minimum tolerance (1e-10)

# IMPORT PACKAGES
from numpy import arange, hstack, vstack, savetxt
#from sympy import N
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
species = "Clavigralla tomentosicollis"
location = "Burkina faso"
period = "Historical"

# USER: Save data to CSV file?
save_data = True


# SELECT INSECT SPECIES
spData = data[data["Species"] == species + " " + location]
#spData = data[data["Species"] == "Clavigralla tomentosicollis Burkina Faso"]
#spData = data[data["Species"] == "Apolygus lucorum China Dafeng"]
#spData = data[data["Species"] == "Adelphocoris suturalis China Dafeng"]
#spData = data[data["Species"] == "Apolygus lucorum China Langfang"]
#spData = data[data["Species"] == "Adelphocoris suturalis China Xinxiang"]
#spData = data[data["Species"] == "Macrosiphum euphorbiae Brazil"]
#spData = data[data["Species"] == "Aulacorthum solani Brazil"]
#spData = data[data["Species"] == "Uroleucon ambrosiae Brazil"]
#spData = data[data["Species"] == "Lygus lineolaris Mississippi"]
#spData = data[data["Species"] == "Pilophorus typicus Japan"]
#spData = data[data["Species"] == "Macrolophus pygmaeus on Myzus persicae Greece"]
#spData = data[data["Species"] == "Macrolophus pygmaeus on Trialeurodes vaporariorum Greece"]


# SELECT LOCATION
#temp_data = temp_data[temp_data["Species"] == species + " " + location]
temp_data = temp_data[temp_data["Species"] == species + " Nigeria"]


# DEFINE MODEL PARAMETERS
# Time parameters
yr = 365 # days in year
init_years = 10 # how many years to use for model initiation
max_years = init_years+80 # how long to run simulations
tstep = 1 # time step = 1 day
CC_years = max_years # how long before climate change "equilibrates"

# Initial abundances
initJ = 1.
initA = 1.
A0 = 0.01 # initial adult density for calculating per capita population growth rate

# Temperature parameters
if period=="Historical":
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
Tmin = spData["Tmin"].values[0] # minimum developmental temperature
# mortality
dJTR = spData["dJTR"].values[0]
AdJ = spData["AdJ"].values[0]
dATR = spData["dATR"].values[0]
AdA = spData["AdA"].values[0]
# competition
qTopt = spData["qTopt"].values[0]
Toptq = spData["Toptq"].values[0]
sq = spData["sq"].values[0]
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

# Life history functions
# fecundity
def b(x):
    #return R(x) * bTopt * exp(-(T(x)-Toptb)**2/2/sb**2)
    return bTopt * exp(-(T(x)-Toptb)**2/2/sb**2)

# maturation rates
def mJ(x):
    return conditional(mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + skew * (exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x))))), 1e-5,
                       1e-5, mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + skew * (exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x)))))) # If mJ(T) < jitcdde min tolerance, then mJ(T) = 0

# mortality rates
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
def M(x):
    return conditional(T(x), Tmin, 0, 1) # if temperature < developmental min, Tmin, then development M = 0; otherwise, M = 1



# DDE MODEL
# Define state variables
J,A,S,τ,r = [y(i) for i in range(5)]


# MODEL
f = {
    J: M(t)*b(t)*A*Allee(A)*exp(-q(t)*A) - M(t)*b(t-τ)*y(1,t-τ)*Allee(y(1,t-τ))*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dJ(t)*J, # juveniles
    
    A: M(t)*b(t-τ)*y(1,t-τ)*Allee(y(1,t-τ))*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dA(t)*A, # Adults
    
    S: S*(mJ(t)/mJ(t-τ)*dJ(t-τ) - dJ(t)), # Through-stage survivorship
    
    τ: 1 - mJ(t)/mJ(t-τ), # Developmental time-delay
    
    r: M(t)*b(t-τ)*A0*Allee(A0)*exp(-q(t-τ)*A0)*mJ(t)/mJ(t-τ)*S - dA(t)*A0 # Low density population growth rate
    }


# RUN DDE SOLVER
# Time and initial conditions
times = arange(0, max_years*yr, tstep)
init = [ initJ, initA, exp(-dJ(-1e-3)/mTR), 1./mTR, 0 ]


# DDE solver
DDE = jitcdde(f, max_delay=1e5, verbose=False)

DDE.constant_past(init)
DDE.compile_C(simplify=False, do_cse=False, verbose=True)
DDE.adjust_diff()


# array containing time and state variables
data = vstack([ hstack([time, DDE.integrate(time)]) for time in times ])
data[:,5] = data[:,5]/data[:,0] # r column from DDE.integrate is actually r*t, so divide by t
data[0,5] = 0 # reset initial r to 0

# SAVE DATA
if save_data == True:
    filename = period + ' time series ' + spData["Species"].values[0] + '.csv'
    savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau,r", comments='')
    #print(DDE.integrate(max_years*yr-180)[:2])
    #print(DDE.integrate(max_years*yr)[:2])


# PLOT
fig,ax = subplots()
ax.plot(data[:,0], data[:,1], label='J')
ax.plot(data[:,0], data[:,2], label='A')
#ax.plot(data[:,0], data[:,3], label='S')
#ax.plot(data[:,0], data[:,4], label='τ')
ax.legend(loc='best')
xlabel("time (days)")
ylabel("population density")
yscale("linear")
xlim((max_years-max_years)*yr,max_years*yr)
ylim(0,100)