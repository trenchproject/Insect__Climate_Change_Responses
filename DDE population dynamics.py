# DELAY DIFFERENTIAL EQUATION MODEL FOR PREDICTING INSECT POPULATION DYNAMICS UNDER SEASONAL TEMPERATURE VARIATION AND CLIMATE CHANGE


# NOTE: If code yields error: "Restarting kernel...", try increasing max_delay in the DDE solver section
# NOTE: If code yields error: "Unsuccessful Integration: Could not integrate with the given tolerance parameters",
#       one of the life history traits is below the minimum tolerance (1e-10)
# NOTE: if code yields error: "CompileError: command 'gcc' failed with exit status 1", one of the parameters is not assigned a value
# NOTE: if code yields error: "IndexError: index 0 is out of bounds for axis 0 with size 0", the location is likely wrong
# NOTE: if code yields warning: "ufunc 'isfinite' not supported for the input types" or "ValueError: cannot convert float NaN to integer",
#       population densities are very high, but this is not a problem for the density-independent model
#       if densities are not plotted, check the output .csv file and remove any rows with 0's after largest density

 
# IMPORT PACKAGES
# Install and import jitcdde from https://github.com/neurophysik/jitcdde
#import sys
#import subprocess
#subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'jitcdde']) # check to see if jitcdde has been installed and if not, install it
from jitcdde import jitcdde, y, t

# Import other necessary packages
from numpy import arange, hstack, vstack, savetxt, isnan
from symengine import exp, pi, cos
from matplotlib.pyplot import subplots, xlabel, ylabel, xlim, ylim, yscale
from pandas import read_csv
from jitcxde_common import conditional
import os


# USER: Enter species, location, and time period (i.e., "Historical" or "Future") or set all_sp to TRUE
species = "Clavigralla tomentosicollis"
location = "Benin"
period = "Historical"
#period = "Future"
all_sp = False

# USER: Save population dynamics data to CSV file?
save_data = True

# USER: Use minimum temperature threshold (i.e., incorporate overwintering)?
overwinter = True

# USER: Include competition (i.e., density-dependent population growth)?
comp = True

# USER: Is model fit to census data?
census = True
    

# DEFINE MODEL PARAMETERS
# Time parameters
yr = 365 # days in year
if comp == True:
    start_yr = 0 # how many years into climate change period to start population dynamics model
    num_yrs = 75 # how long to run model
else: # NOTE: density-independent model is run between years 65 and 75 of climate change because otherwise densities become too large due to unbounded population growth
    start_yr = 65
    num_yrs = 10
if census == True:
    start_yr = 0
    num_yrs = 10
    if species == "Apolygus lucorum":
        active_date = 137 # NOTE: for Apolygus lucorum during the census, insects only become active once the fields were planted on in mid-May (day = 137)
tstep = 1 # time step = 1 day
    
# Initial abundances
initJ = 100.
initA = 10.


# INPUT TEMPERATURE RESPONSE PARAMETERS AND TEMPERATURE PARAMETERS
# Set path to temperature response and habitat temperature files
path = os.getcwd()
for root, dirs, files in os.walk(os.getcwd()):
        if "Johnson_Insect_Responses" in dirs:
            path = os.path.join(root, "Johnson_Insect_Responses")
TR_params = read_csv(path + "/Model parameters/Temperature response parameters.csv")
T_params = read_csv(path + "/Model parameters/Habitat temperature parameters.csv")



# SELECT SPECIES
if all_sp == True:
    sp = 0
else:
    sp = TR_params[TR_params["Species"] == species + " " + location].index[0]


# REPEAT CODE FOR ALL SPECIES
while(True):
    
    # SELECT SPECIES AND ASSIGN PARAMETERS
    sp_params = TR_params.iloc[sp]
    t_params = T_params.iloc[sp]
    
    # Temperature parameters (see "Habitat temperature parameters.csv")
    if period == "Historical":
        meanT = t_params["meanT.h"]
        amplT = t_params["amplT.h"] 
        shiftT = t_params["shiftT.h"]
        delta_mean = t_params["delta_mean.h"]
        delta_ampl = t_params["delta_ampl.h"]
    else:
        meanT = t_params["meanT.f"]
        amplT = t_params["amplT.f"] 
        shiftT = t_params["shiftT.f"]
        delta_mean = t_params["delta_mean.f"]
        delta_ampl = t_params["delta_ampl.f"]
    
    # Life history and competitive traits (see "Temperature response parameters.csv")
    # fecundity
    bTopt = sp_params["bTopt"]
    Toptb = sp_params["Toptb"]
    sb = sp_params["sb"]
    # maturation
    gMax = sp_params["gMax"]
    gTR = sp_params["gTR"]
    TR = sp_params["TR"]
    Ag = sp_params["Ag"]
    AL = sp_params["AL"]
    TL = sp_params["TL"]
    Tmin = sp_params["Tmin"]
    Toptg = sp_params["Toptg"]
    Tmaxg = sp_params["Tmaxg"]
    # mortality
    dJTR = sp_params["dJTR"]
    AdJ = sp_params["AdJ"]
    dATR = sp_params["dATR"]
    AdA = sp_params["AdA"]
    # competition
    if comp == True:
        qTopt = sp_params["qTopt"]
    else:
        qTopt = 0
    Toptq = Toptb
    sq = sb
    
    
    # FUNCTIONS
    # Seasonal temperature variation (Eq. 5)
    def T(x):
            return (meanT + delta_mean*(x+start_yr*yr)) - (amplT + delta_ampl*(x+start_yr*yr)) * cos(2*pi*(x + shiftT)/yr)
    
    # Life history functions
    # per capita birth rate (Eq. 2a; see Eq. S3a for overwintering)
    def b(x):
        return conditional(bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2)), 1e-5, 1e-5, bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2))) # If b(T) < jitcdde min tolerance, then b(T) = 0
    
    # stage-specific per capita mortality rates (Eq. 2b; see Eq. S3b for overwintering)
    def dJ(x):
        return dJTR * exp(AdJ * (1/TR - 1/T(x)))
    def dA(x):
        return conditional(T(x), Tmin, 0.1 + dATR * exp(AdA * (1/TR - 1/T(x))), dATR * exp(AdA * (1/TR - 1/T(x)))) # if habitat temperature < developmental min (Tmin), then dA = dA + 0.1; otherwise, use dA(T(x))
   
    # development rate (Eq. 2c)
    def g(x):
        return conditional(gTR * T(x)/TR * exp(Ag * (1/TR - 1/T(x))) / (1 + exp(AL*(1/TL-1/T(x)))), 1e-5, 1e-5, # If g(T) < jitcdde min tolerance, then g(T) = 1e-5
                               conditional(T(x), Toptg, gTR * T(x)/TR * exp(Ag * (1/TR - 1/T(x))) / (1 + exp(AL*(1/TL-1/T(x)))), # If habitat temperature < developmental optima (Toptg), then use monotonic g(T)
                               conditional(T(x), Tmaxg, gMax, 1e-5)))  # If habitat temperature < developmental maximum (Tmaxg), then g(T) = gMax; otherwise, g(T) = 1e-5
        
    # density-dependence due to competition
    def q(x):
        return qTopt * exp(-(T(x)-Toptq)**2/(2*sq**2))
    
    # Overwintering
    if overwinter == False:
        def O(x):
            return 1
    else:
        if census == False or species == "Clavigralla tomentosicollis":
            def O(x):
                return conditional(T(x), Tmin, 0, 1) # if habitat temperature < developmental min (Tmin), then O = 0 (i.e., insect is overwintering); otherwise, O = 1 (i.e., insect is active)
        else:
            def O(x):
                return conditional(T(x), T(active_date), 0, 1) # if habitat temperature < temperature in mid-May, then O = 0 (i.e., insect is overwintering); otherwise, O = 1 (i.e., insect is active)

    
    # DDE MODEL
    # Define state variables
    J,A,S,τ = [y(i) for i in range(4)]
    
    # DDE model
    f = {
        J: O(t)*b(t)*A*exp(-q(t)*A) - O(t-τ)*b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*g(t)/g(t-τ)*S - dJ(t)*J, # juvenile density (Eq. 3a incorporating Eq. 4a); NOTE: python names this variable "y(0)"
        
        A: O(t-τ)*b(t-τ)*y(1,t-τ)*exp(-q(t-τ)*y(1,t-τ))*g(t)/g(t-τ)*S - dA(t)*A, # Adult density (Eq. 3b incorporating Eq. 4a); NOTE: python names this variable "y(1)"
        
        S: S*(g(t)/g(t-τ)*dJ(t-τ) - dJ(t)), # Juvenile through-stage survival (Eq. 4c); NOTE: python names this variable "y(2)"
        
        τ: 1 - g(t)/g(t-τ), # Developmental time-delay (Eq. 4b); NOTE: python names this variable "y(3)"
        }
    
    
    # RUN DDE SOLVER
    # Time and initial conditions
    times = arange(0, num_yrs*yr, tstep)
    init = [ initJ, initA, exp(-dJ(-1e-3)/g(-1e-3)), 1./g(-1e-3)]
    
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
        if census == False:
            if comp == True:
                filename = 'Time series data new/' + period + ' time series ' + sp_params["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            else:
                filename = 'Time series data DI/' + period + ' time series ' + sp_params["Species"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
        else:
            filename = 'Time series data Census new/Census time series ' + sp_params["Species"] + '.csv'
            savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
    
    
    # PLOT
    fig,ax = subplots()
    ax.plot(data[:,0], data[:,1], label='J')
    ax.plot(data[:,0], data[:,2], label='A')
    ax.plot(data[:,0], vstack([100*O(i) for i in arange(0,num_yrs*yr,1) ]), label='overwinter')
    ax.plot(data[:,0], vstack([T(i) - 273 for i in arange(0,num_yrs*yr,1) ]), label='T')
    #ax.plot(data[:,0], data[:,3], label='S')
    #ax.plot(data[:,0], data[:,4], label='τ')
    ax.legend(loc='best')
    xlabel("time (days)")
    ylabel("population density")
    yscale("linear")
    xlim((num_yrs-2)*yr,(num_yrs-0)*yr)
    ylim(0,200)
    
    
    # END LOOP WHEN MODEL IS RUN FOR ALL SPECIES
    if all_sp == True:
        sp = sp + 1
    else:
        sp = TR_params.shape[0] # if all_sp == False, then set "sp" to a value that will cause the loop to break
    if sp > TR_params.shape[0] - 1:
        break
