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


# USER: Save data to CSV file?
save_data = True

# USER: Use minimum temperature threshold?
minT = True

# USER: Incorporate diurnal temperature fluctuations?
daily = False


# INPUT TEMPERATURE RESPONSE PARAMETERS AND TEMPERATURE PARAMETERS
S_data = read_csv("Temperature response parameters.csv")
if daily == True:
    T_data = read_csv("Temperature parameters.csv")
else:
    T_data = read_csv("Temperature parameters Tave.csv")


# REPEAT CODE FOR ALL SPECIES
species = 15
while(True):
    
    # SELECT SPECIES
    spData = S_data.iloc[species]
    temp_data = T_data.iloc[species]
        
    
    # DEFINE MODEL PARAMETERS
    # Time parameters
    yr = 365 # days in year
    start_date = 0 # day on which to start model
    init_years = 0 # how many years into climate change to start model
    max_years = init_years + 125 # how long to run simulations
    tstep = 1 # time step = 1 day
    
    # Initial abundances
    initJ = 100.
    initA = 100.
    
    # Temperature parameters
    meanT = temp_data["meanT.h"]
    amplT = temp_data["amplT.h"] 
    shiftT = temp_data["shiftT.h"]
    delta_mean = temp_data["delta_mean.h"]
    delta_ampl = temp_data["delta_ampl.h"]
    amplD = temp_data["amplD.h"]
    
    # Increase mean temperature by 0.1 until extinct
    delta_mean = 0.1/yr
    
    # Life history and competitive traits
    # fecundity
    bTopt = spData["bTopt"]
    Toptb = spData["Toptb"]
    sb = spData["sb"]
    # maturation
    mTR = spData["mTR"]
    TR = spData["TR"]
    AmJ = spData["AmJ"]
    skew = spData["skew"]
    AL = spData["AL"]
    TL = spData["TL"]
    AH = spData["AH"]
    TH = spData["TH"]
    Tmin = spData["Tmin"] # minimum developmental temperature
    # mortality
    dJTR = spData["dJTR"]
    AdJ = spData["AdJ"]
    dATR = spData["dATR"]
    AdA = spData["AdA"]
    # competition
    qTopt = spData["qTopt"]
    Toptq = Toptb #spData["Toptq"]
    sq = sb #spData["sq"]
    #Aq = spData["Aq"]
    #Tmax = spData["Tmax"]
    #qTopt = qTR*exp(Aq*(1/TR - 1/Tmax))
    
    
    # FUNCTIONS
    # Seasonal temperature variation (K) over time
    def T(x):
            return conditional(x, start_date, meanT, # during "pre-history", habitat temperature is constant at its mean
                               (meanT + delta_mean*(x+init_years*yr)) - (amplT + delta_ampl*(x+init_years*yr)) * cos(2*pi*(x + shiftT)/yr)  - amplD * cos(2*pi*x)) # temperature regime during climate change
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
        if daily == True:
            filename = 'Time series data diurnal Ext/Time series ' + spData["Species"] + '.csv'
            savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
        else:
            filename = 'Time series data Ext/Time series ' + spData["Species"] + '.csv'
            savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')


    # END LOOP WHEN MODEL IS RUN FOR ALL SPECIES
    species = species + 1
    if species > S_data.shape[0] - 1:
        break


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
xlim((max_years-max_years)*yr,(max_years-0)*yr)
ylim(0,200)
