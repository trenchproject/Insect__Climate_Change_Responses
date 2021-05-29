# Delay differential equation model for predicting insect population dynamics
# under seasonal temperature variation and climate change

# IMPORT PACKAGES
from numpy import array, arange, hstack, vstack, savetxt, pi
from jitcdde import jitcdde, y, t
from symengine import exp, sin
from matplotlib.pyplot import subplots, xlabel, ylabel, yscale, xlim, ylim
from pandas import read_csv


# INPUT TEMPERATURE RESPONSE DATA
tempData = read_csv("Temperature response data.csv")


# SELECT INSECT SPECIES
spData = tempData[tempData["Species"] == "Clavigralla shadabi"]
#spData = tempData[tempData["Species"] == "Clavigralla tomentosicollis"]   


# DEFINE MODEL PARAMETERS
# Habitat temperature and climate change parameters
meanT = spData["meanT"].values[0]
amplT = spData["amplT"].values[0] 
shiftT = spData["shiftT"].values[0]
delta_mean = 0.
delta_ampl = 0.

# Time parameters
yr = 365.0 # days in year
max_years = 100. # how long to run simulations
tstep = 1. # time step = 1 day
delta_years = 100 # how long before climate change "equilibrates"

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
# mortality
dJTR = spData["dJTR"].values[0]
AdJ = spData["AdJ"].values[0]
dATR = spData["dATR"].values[0]
AdA = spData["AdA"].values[0]
# competition
qTR = spData["qTR"].values[0]
Toptq = spData["Toptq"].values[0]
sq = spData["sq"].values[0]
Aq = spData["Aq"].values[0]
Tmax = spData["Tmax"].values[0]
qTemp = spData["qTemp"].values[0]
qTopt = qTR*exp(Aq*(1./TR - 1./Tmax))


# FUNCTIONS
# Seasonal temperature variation (K) over time
def T(t):
    # climate change: climate will change by delta_mean and delta_ampl over delta_years
    m_mean = (delta_mean/delta_years)*(1/yr)
    m_ampl = (delta_ampl/delta_years)*(1/yr)
    if t < 0:
        return meanT # during model initiation, habitat temperature is constant at its mean
    elif t < (delta_years*yr):
        return (meanT + m_mean*t) + (amplT + m_ampl*t) * sin(2*pi*(t + shiftT)/yr) # temperature variation during climate change
    else:
        return (meanT + delta_mean) + (amplT + delta_ampl) * sin(2*pi*(t + shiftT)/yr) # temperature varation after climate change "equilibriates"


# Life history functions
# fecundity
def b(t):
    return bTopt * exp(-(T(t)-Toptb)**2/2./sb**2)

# maturation rates
def mJ(t):
    return mTR * T(t)/TR * exp(AmJ * (1./TR - 1./T(t))) / (1. + skew * (exp(AL*(1./TL-1./T(t)))+exp(AH*(1./TH-1./T(t)))))

# mortality rates
def dJ(t):
    return dJTR * exp(AdJ * (1./TR - 1./T(t)))
def dA(t):
    return dATR * exp(AdA * (1./TR - 1./T(t)))

# density-dependence due to competition
def q(t):
    return (1-qTemp) * qTR + qTemp * qTopt * exp(-(T(t)-Toptq)**2/2./sq**2) # q is temperature dependent if qTemp = 1 or constant if qTemp = 0



# DDE MODEL
# Define state variables
J,A,S,τ = [y(i) for i in range(4)]
labels = ["J","A","S","τ"]


# Model
f = {
    J: b(t)*A*exp(-q(t)*A) - b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dJ(t)*J,
    
    A: b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dA(t)*A,
    
    S: S*(mJ(t)/mJ(t-τ)*dJ(t-τ) - dJ(t)),
    
    τ: 1. -  mJ(t)/mJ(t-τ)
    }



# RUN DDE SOLVER
# Time and initial conditions
times = arange(0.0, max_years*yr, tstep)
init = array([ 1., 1., exp(-dJ(-1e-3)/mTR), 1./mTR ])


# DDE solver
DDE = jitcdde(f, max_delay=1e10)
DDE.constant_past(init)
DDE.integrate_blindly(0.01)


# Save data array containing time and state variables
data = vstack([ hstack([time,DDE.integrate(time)]) for time in times ])
savetxt("timeseries.dat", data) #, fmt='%s')



# PLOT
fig,ax = subplots()
ax.plot(data[:,0], data[:,1], label='J')
ax.plot(data[:,0], data[:,2], label='A')
ax.legend(loc='best')
xlabel("time (days)")
ylabel("population density")
yscale("log")
xlim((max_years-1)*yr,max_years*yr)
ylim(0.1,100)