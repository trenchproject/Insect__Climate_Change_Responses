# Delay differential equation model for predicting insect population dynamics
# under seasonal temperature variation and climate change
# NOTE: If code yields error: "Restarting kernel...", try increasing max_delay in the DDE solver section
# NOTE: If code yields error: "Unsuccessful Integration: Could not integrate with the given tolerance parameters",
#       one of the life history traits is below the minimum tolerance (1e-10)

# IMPORT PACKAGES
from numpy import arange, hstack, vstack, savetxt
#from sympy import N
from jitcdde import jitcdde, y, t
from symengine import exp, pi, sin, cos, asin
from matplotlib.pyplot import subplots, xlabel, ylabel, xlim, ylim, yscale, plot, show
from pandas import read_csv
from jitcxde_common import conditional


# INPUT TEMPERATURE RESPONSE DATA
tempData = read_csv("Temperature response parameters.csv")


# SELECT INSECT SPECIES
#spData = tempData[tempData["Species"] == "Clavigralla shadabi Benin"]
#spData = tempData[tempData["Species"] == "Clavigralla tomentosicollis Benin"]   
#spData = tempData[tempData["Species"] == "Clavigralla tomentosicollis Nigeria B"]
#spData = tempData[tempData["Species"] == "Clavigralla tomentosicollis Burkina Faso"]
spData = tempData[tempData["Species"] == "Apolygus lucorum China Dafeng"]
#spData = tempData[tempData["Species"] == "Adelphocoris suturalis China Dafeng"]
#spData = tempData[tempData["Species"] == "Apolygus lucorum China Langfang"]
#spData = tempData[tempData["Species"] == "Adelphocoris suturalis China Xinxiang"]
#spData = tempData[tempData["Species"] == "Macrosiphum euphorbiae Brazil"]
#spData = tempData[tempData["Species"] == "Aulacorthum solani Brazil"]
#spData = tempData[tempData["Species"] == "Uroleucon ambrosiae Brazil"]
#spData = tempData[tempData["Species"] == "Lygus lineolaris Mississippi"]
#spData = tempData[tempData["Species"] == "Pilophorus typicus Japan"]
#spData = tempData[tempData["Species"] == "Macrolophus pygmaeus on Myzus persicae Greece"]
#spData = tempData[tempData["Species"] == "Macrolophus pygmaeus on Trialeurodes vaporariorum Greece"]

# DEFINE MODEL PARAMETERS
# Time parameters
yr = 360 # days in year
init_years = 10 # how many years to use for model initiation
max_years = init_years+140 # how long to run simulations
tstep = 1 # time step = 1 day
CC_years = max_years # how long before climate change "equilibrates"

# Initial abundances
initJ = 1.
initA = 1.

# Habitat temperature and climate change parameters
meanT = spData["meanT"].values[0]
amplT = spData["amplT"].values[0] 
shiftT = spData["shiftT"].values[0]
delta_mean = spData["delta_mean"].values[0]
delta_ampl = spData["delta_ampl"].values[0]
#delta_mean = 2.5/(100*yr)
#delta_ampl = 2.5/(100*yr)

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
qTR = spData["qTR"].values[0]
Toptq = spData["Toptq"].values[0]
sq = spData["sq"].values[0]
Aq = spData["Aq"].values[0]
Tmax = spData["Tmax"].values[0]
qTemp = spData["qTemp"].values[0]
qTopt = qTR#*exp(Aq*(1/TR - 1/Tmax))
'''
# Resource variation
Res = spData["Res"].values[0] # Res = 1 (0) if there is (not) resource variation
Rstart = spData["Rstart"].values[0] # start date of resource availablity
Rend = spData["Rend"].values[0] # end of resource availability
if Rend > Rstart:               # proportion of year when resource is available
    Rlength = (Rend - Rstart)/yr
else:
     Rlength = (yr - Rstart + Rend)/yr
Rshift = cos(pi*Rlength) # shift used to model resource availability via sine wave
'''


# FUNCTIONS

# Seasonal temperature variation (K) over time
def T(x):
    return conditional(x, 0, meanT, # during "pre-history" (t<0), habitat temperature is constant at its mean
                       conditional(x, init_years*yr, meanT + amplT * sin(2*pi*(x + shiftT)/yr), # during model initiation, delta_mean = 0 and delta_ampl = 0
                                   conditional(x, CC_years*yr, (meanT + delta_mean*(x-init_years*yr)) + (amplT + delta_ampl*(x-init_years*yr)) * sin(2*pi*((x-init_years*yr) + shiftT)/yr), # temperature regime during climate change
                                               (meanT + delta_mean*CC_years*yr) + (amplT + delta_ampl*CC_years*yr) * sin(2*pi*((x-init_years*yr) + shiftT)/yr)))) # temperature regime after climate change "equilibriates"
'''
# Plot temperature function
xvals = arange(0,10*yr,1)
yvals = vstack([(meanT + delta_mean*(i-init_years*yr)) + (amplT + delta_ampl*(i-init_years*yr)) * sin(2*pi*(i-init_years*yr + shiftT)/yr) for i in xvals ])
plot(xvals,yvals)
show()
'''
'''
# Seasonal temperature variation (K) over time based on Fourier analysis of historical data
def T(x):
    return conditional(x, 0, meanT, # during "pre-history" (t<0), habitat temperature is constant at its mean
                       conditional(x, init_years*yr, 299.08 + 1.84*cos(2*pi*30/(30*yr)*(x-15) - 1.86) + 1.51*cos(2*pi*60/(30*yr)*(x-15) - 3.09) + 0.359*cos(2*pi*90/(30*yr)*(x-15) - 2.57), # during model initiation, no change in mean temperature
                                   conditional(x, CC_years*yr, 299.08 + 0*0.000127*x + delta_mean*(x-init_years*yr) + (1.84*cos(2*pi*30/(30*yr)*(x-init_years*yr-15) - 1.86) + 1.51*cos(2*pi*60/(30*yr)*(x-init_years*yr-15) - 3.09) + 0.359*cos(2*pi*90/(30*yr)*(x-init_years*yr-15) - 2.57)), # temperature regime during climate change
                                               299.08 + 0.000127*CC_years*yr + delta_mean*CC_years*yr + 1.84*cos(2*pi*30/(30*yr)*(x-init_years*yr-15) - 1.86) + 1.51*cos(2*pi*60/(30*yr)*(x-init_years*yr-15) - 3.09) + 0.359*cos(2*pi*90/(30*yr)*(x-init_years*yr-15) - 2.57)))) # temperature regime after climate change "equilibriates"

# Plot temperature function
xvals = arange(0,1*yr,1)
yvals = vstack([299.08 + (1.84*cos(2*pi*30/(30*yr)*(i-init_years*yr-15) - 1.86) + 1.51*cos(2*pi*60/(30*yr)*(i-init_years*yr-15) - 3.09) + 0.359*cos(2*pi*90/(30*yr)*(i-init_years*yr-15) - 2.57)) for i in xvals ])
plot(xvals,yvals)
show()
'''
'''
# Seasonal variation in resource quality (Res = 1: resource available, Res = 0: resource unavailable)
def R(x):
    dummy = 1/(1-Rshift)*(-Rshift + sin(2*pi*(x-Rstart)/yr + asin(Rshift))) # 'dummy' sine wave describing resource availability
    if Res == 1: # incorporate resource variation in model
        return 0.5*(abs(dummy)/dummy + 1) # set negative phases of sine wave to zero and positive phases to 1
        #return 0.5*(abs(dummy) + dummy) # set negative phases of sine wave to zero
    else:
        return 1 # do not incorporate resource variation in model

# Plot resource function
xvals = arange(0,2*yr,1)
yvals = vstack([0.5*(abs(1/(1-Rshift)*(-Rshift + sin(2*pi*(i-Rstart)/yr + asin(Rshift)))) + 1/(1-Rshift)*(-Rshift + sin(2*pi*(i-Rstart)/yr + asin(Rshift)))) for i in xvals ])
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
    return mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + skew * (exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x)))))

# mortality rates
def dJ(x):
    return dJTR * exp(AdJ * (1/TR - 1/T(x)))
def dA(x):
    #dummy = R(x)/(abs(R(x)-1)+1) # 'dummy' function that is 1 when R(x) > 0 and 0 when R(x) = 0
    #return 0.02*(1 - dummy) + dummy * dATR * exp(AdA * (1/TR - 1/T(x)))
    return dATR * exp(AdA * (1/TR - 1/T(x)))

# density-dependence due to competition
def q(x):
    return (1-qTemp) * qTopt + qTemp * qTopt * exp(-(T(x)-Toptq)**2/2/sq**2) # q is temperature dependent if qTemp = 1 or constant if qTemp = 0

# Allee effect
A_thr = 0.00*qTopt # Allee threshold
def Allee(x):
    return conditional(x, A_thr, 0, 1) # if A <= A_thr then Allee = 0, otherwise Allee = 1 


# DDE MODEL
# Define state variables
J,A,S,τ = [y(i) for i in range(4)]


# MODEL
f = {
    J: b(t)*A*Allee(A)*exp(-q(t)*A) - b(t-τ)*y(1,t-τ)*Allee(y(1,t-τ))*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dJ(t)*J, # juvenile density
    
    A: b(t-τ)*y(1,t-τ)*Allee(y(1,t-τ))*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dA(t)*A, # Adult density
    
    S: S*(mJ(t)/mJ(t-τ)*dJ(t-τ) - dJ(t)), # Through-stage survivorship
    
    τ: 1 -  mJ(t)/mJ(t-τ) # Developmental time-delay
    }


# RUN DDE SOLVER
# Time and initial conditions
times = arange(0, max_years*yr, tstep)
init = [ initJ, initA, exp(-dJ(-1e-3)/mTR), 1./mTR ]


# DDE solver
DDE = jitcdde(f, max_delay=1e5, verbose=False)

DDE.constant_past(init)
DDE.compile_C(simplify=False, do_cse=False, verbose=True)
DDE.adjust_diff()


# Save data array containing time and state variables
data = vstack([ hstack([time, DDE.integrate(time)]) for time in times ])
filename = 'Time series ' + spData["Species"].values[0] + '.csv'  
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
yscale("linear")
xlim((max_years-max_years)*yr,max_years*yr)
ylim(0,100)