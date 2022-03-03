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
from matplotlib.pyplot import subplots, xlabel, ylabel, xlim, ylim, yscale, plot
from pandas import read_csv
from jitcxde_common import conditional
import os


# SET WORKING DIRECTORY TO SAME AS THE PYTHON CODE
cwd = os.getcwd()
if cwd != '/Users/johnson/Documents/Christopher/GitHub/Johnson_Insect_Responses':
    os.chdir('/Users/johnson/Documents/Christopher/GitHub/Johnson_Insect_Responses')


# USER: Enter species, location, and time period
species = "Brevicoryne brassicae"
location = "US Columbia"
period = "Historical"
#period = "Future"

# USER: Run model for all species?
all_sp = True

# USER: Save data to CSV file?
save_data = True

# USER: Use minimum temperature threshold?
minT = True

# USER: Use left-skewed function for development?
left_skew = False # if False, development plateaus between Topt and Tmax before going to zero above Tmax

# USER: Include competition (i.e., density-dependent population growth)?
comp = True

# USER: Incorporate diurnal temperature fluctuations?
daily = False

# USER: Is model fit to census data?
census = False
    

# DEFINE MODEL PARAMETERS
# Time parameters
yr = 365 # days in year
start_date = 0 # day on which to start model
if comp == True:
    init_years = 0 # how many years into climate change to start model
    max_years = start_date + 75 # how long to run simulations
else:
    init_years = 65
    max_years = start_date + 10 # NOTE: no model initialization (start_date = 0) because densities become too large before model begins due to unbounded population growth
if census == True:
    init_years = 65 # how many years into climate change to start model
    max_years = start_date + 10 # how long to run simulations
    census_start = 137 # for Apolygus lucorum, start date is May 17 (day = 137) because the field was planted in mid-May
    rainy_season = 306 # for Clavigralla tomentosicollis, rainy season start date is Nov 1 (day = 306)
tstep = 1 # time step = 1 day
CC_years = max_years # how long before climate change "equilibrates"
    
# Initial abundances
initJ = 100.
initA = 100.


# INPUT TEMPERATURE RESPONSE PARAMETERS AND TEMPERATURE PARAMETERS
Sdata = read_csv("Temperature response parameters.csv")
if daily == True:
    Tdata = read_csv("Temperature parameters.csv")
else:
    Tdata = read_csv("Temperature parameters Tave.csv")


# SELECT SPECIES
if all_sp == True:
    sp = 0
else:
    sp = Sdata[Sdata["Species"] == species + " " + location].index[0]


# REPEAT CODE FOR ALL SPECIES
while(True):
    
    # SELECT SPECIES
    spData = Sdata.iloc[sp]
    tempData = Tdata.iloc[sp]

    
    # Temperature parameters
    if period == "Historical":
        meanT = tempData["meanT.h"]
        amplT = tempData["amplT.h"] 
        shiftT = tempData["shiftT.h"]
        delta_mean = tempData["delta_mean.h"]
        delta_ampl = tempData["delta_ampl.h"]
        amplD = tempData["amplD.h"]
    else:
        meanT = tempData["meanT.f"]
        amplT = tempData["amplT.f"] 
        shiftT = tempData["shiftT.f"]
        delta_mean = tempData["delta_mean.f"]
        delta_ampl = tempData["delta_ampl.f"]
        amplD = tempData["amplD.f"]
    if daily == False:
        amplD = 0  
    
    # Life history and competitive traits
    # fecundity
    bTopt = spData["bTopt"]
    Toptb = spData["Toptb"]
    sb = spData["sb"]
    # maturation
    mTR = spData["mTR"]
    TR = spData["TR"]
    AmJ = spData["AmJ"]
    AL = spData["AL"]
    TL = spData["TL"]
    AH = spData["AH"]
    TH = spData["TH"]
    Tmin = spData["Tmin"] # minimum developmental temperature
    Topt = spData["Topt"] # optimum developmental temperature
    Tmax = spData["Tmax"] # maximum developmental temperature
    # mortality
    dJTR = spData["dJTR"]
    AdJ = spData["AdJ"]
    dATR = spData["dATR"]
    AdA = spData["AdA"]
    # competition
    if comp == True:
        qTopt = spData["qTopt"]
    else:
        qTopt = 0
    Toptq = Toptb #spData["Toptq"]
    sq = sb #spData["sq"]
    #Aq = spData["Aq"]
    #Tmax = spData["Tmax"]
    #qTopt = qTR*exp(Aq*(1/TR - 1/Tmax))
    
    
    # FUNCTIONS
    # Seasonal temperature variation (K) over time
    def T(x):
            return conditional(x, start_date*yr, meanT, # during "pre-history", habitat temperature is constant at its mean
                               conditional(x, CC_years*yr, (meanT + delta_mean*(x+init_years*yr)) - (amplT + delta_ampl*(x+init_years*yr)) * cos(2*pi*(x + shiftT)/yr)  - amplD * cos(2*pi*x), # temperature regime during climate change
                                           (meanT + delta_mean*CC_years*yr) - (amplT + delta_ampl*CC_years*yr) * cos(2*pi*(x + shiftT)/yr)  - amplD * cos(2*pi*x))) # temperature regime after climate change "equilibriates"
    '''
    # Plot temperature function
    xvals = arange(0,1*yr,0.1)
    yvals = vstack([(meanT + delta_mean*i) - (amplT + delta_ampl*i) * cos(2*pi*(i + shiftT)/yr)  - amplD * cos(2*pi*i) for i in xvals ])
    plot(xvals,yvals)
    '''
    
    # Life history functions
    # fecundity
    def b(x):
        return conditional(bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2)), 1e-5, 1e-5, bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2))) # If b(T) < jitcdde min tolerance, then b(T) = 0
    
    # maturation rates
    if left_skew == True:
        def mJ(x):
            return conditional(mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + 0*(exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x))))), 1e-5,
                               1e-5, mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + 0*(exp(AL*(1/TL-1/T(x)))+exp(AH*(1/TH-1/T(x)))))) # If mJ(T) < jitcdde min tolerance, then mJ(T) = 1e-5
    else:
        def mJ(x):
            return conditional(mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + exp(AL*(1/TL-1/T(x)))), 1e-5, 1e-5, # If mJ(T) < jitcdde min tolerance, then mJ(T) = 1e-5
                               conditional(T(x), Topt, mTR * T(x)/TR * exp(AmJ * (1/TR - 1/T(x))) / (1 + exp(AL*(1/TL-1/T(x)))), # If temperature < developmental optima (Topt), then use monotonic mJ(T)
                               conditional(T(x), Tmax, mTR * Topt/TR * exp(AmJ * (1/TR - 1/Topt)) / (1 + exp(AL*(1/TL-1/Topt))), 1e-5)))  # If temperature < developmental maximum (Tmax), then mJ(T) = mJ(Topt); otherwise, mJ(T) = 1e-5
        
    # mortality rates
    def dJ(x):
        return dJTR * exp(AdJ * (1/TR - 1/T(x)))  
    def dA(x):
        return conditional(T(x), Tmin, 0.1 + dATR * exp(AdA * (1/TR - 1/T(x))), dATR * exp(AdA * (1/TR - 1/T(x)))) # if temperature < developmental min (Tmin), then dA = dA + 0.1; otherwise, use dA(T(x))
   
    # density-dependence due to competition
    def q(x):
        return qTopt * exp(-(T(x)-Toptq)**2/(2*sq**2))
    
    
    # Overwintering in temperate habitats and seasonal resource variability in tropical and subtropical habitats
    if minT == False:
        def R(x):
            return 1
        def M(x):
            return 1
    else:
        if census == False:
            def R(x):
                return 1
            def M(x):
                return conditional(T(x), Tmin, 0, 1) # if temperature < developmental min (Tmin) then development M = 0; otherwise, M = 1
        else:
            if species == "Clavigralla tomentosicollis":
                def R(x):
                    return conditional(sin(2*pi*(x-rainy_season)/365), 0, 0, 1)
                def M(x):
                    return 1
            else:
                def R(x):
                    return 1
                if (species == "Apolygus lucorum" or species == "Adelphocoris suturalis") and period == "Historical":
                    def M(x):
                        return conditional(T(x), T(census_start), 0, 1)
                else:
                    def M(x):
                        return conditional(T(x), Tmin, 0, 1)

    
    # DDE MODEL
    # Define state variables
    J,A,S,τ = [y(i) for i in range(4)]
    
    # DDE model
    f = {
        J: conditional(t, start_date, 0, M(t)*R(t)*b(t)*A*exp(-q(t)*A) - M(t-τ)*R(t-τ)*b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dJ(t)*J), # juveniles: y(0)
        
        A: conditional(t, start_date, 0, M(t-τ)*R(t-τ)*b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*mJ(t)/mJ(t-τ)*S - dA(t)*A), # Adults: y(1)
        
        S: conditional(t, start_date, 0, S*(mJ(t)/mJ(t-τ)*dJ(t-τ) - dJ(t))), # Through-stage survivorship: y(2)
        
        τ: conditional(t, start_date, 0, 1 - mJ(t)/mJ(t-τ)), # Developmental time-delay: y(3)
        }
    
    
    # RUN DDE SOLVER
    # Time and initial conditions
    times = arange(0, max_years*yr, tstep)
    init = [ initJ, initA, exp(-dJ(-1e-3)/mJ(-1e-3)), 1./mJ(-1e-3)]
    
    # Run DDE solver
    DDE = jitcdde(f, max_delay=1e5, verbose=False)
    DDE.constant_past(init)
    DDE.compile_C(simplify=False, do_cse=False, verbose=True)
    DDE.adjust_diff()
    
    
    # SAVE DATA
    # array containing time and state variables (can add N(T(time)) to output temperature data, but with significantly longer run-time)
    data = vstack([ hstack([time, DDE.integrate(time)]) for time in times ])
    
    # set values below 1e-5 or NAN to 0
    data[data < 1e-5] = 0
    data[isnan(data)] = 0
    
    # save data to csv 
    if save_data == True:
        if census == False:
            if comp == True and daily == True and left_skew == True:
                filename = 'Time series data Diurnal/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if comp == True and daily == False and left_skew == True:
                filename = 'Time series data Tave Mon/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if comp == False and daily == True and left_skew == True:
                filename = 'Time series data DI Diurnal/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if comp == False and daily == False and left_skew == True:
                filename = 'Time series data DI Tave/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if comp == True and daily == True and left_skew == False:
                filename = 'Time series data Diurnal Dev/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if comp == True and daily == False and left_skew == False:
                filename = 'Time series data Tave Dev/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if comp == False and daily == True and left_skew == False:
                filename = 'Time series data DI Diurnal Dev/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if comp == False and daily == False and left_skew == False:
                filename = 'Time series data DI Tave Dev/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
        # For fits to census data
        else:
            if daily == True and left_skew == True:
                filename = 'Time series data Census/' + period + ' time series ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if daily == False and left_skew == True:
                filename = 'Time series data Census/' + period + ' time series Tave ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if daily == True and left_skew == False:
                filename = 'Time series data Census/' + period + ' time series Dev ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            if daily == False and left_skew == False:
                filename = 'Time series data Census/' + period + ' time series Tave Dev ' + spData["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
    
    
    # PLOT
    fig,ax = subplots()
    ax.plot(data[:,0], data[:,1], label='J')
    ax.plot(data[:,0], data[:,2], label='A')
    ax.plot(data[:,0], vstack([100*M(i) for i in arange(0,max_years*yr,1) ]), label='M')
    ax.plot(data[:,0], vstack([T(i) - 273 for i in arange(0,max_years*yr,1) ]), label='T')
    #ax.plot(data[:,0], data[:,3], label='S')
    #ax.plot(data[:,0], data[:,4], label='τ')
    ax.legend(loc='best')
    xlabel("time (days)")
    ylabel("population density")
    yscale("linear")
    xlim((max_years-2)*yr,(max_years-0)*yr)
    ylim(0,200)
    
    
    # END LOOP WHEN MODEL IS RUN FOR ALL SPECIES
    if all_sp == True:
        sp = sp + 1
    else:
        sp = Sdata.shape[0]
    if sp > Sdata.shape[0] - 1:
        break
