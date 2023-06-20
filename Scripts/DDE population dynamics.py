# DELAY DIFFERENTIAL EQUATION MODEL FOR PREDICTING INSECT POPULATION DYNAMICS UNDER SEASONAL TEMPERATURE VARIATION AND CLIMATE CHANGE


# NOTE: If code yields error: "Restarting kernel...", try increasing max_delay in the DDE solver section
# NOTE: If code yields error: "Unsuccessful Integration: Could not integrate with the given tolerance parameters",
#       one of the life history traits is below the minimum tolerance (1e-10)
# NOTE: if code yields error: "CompileError: command 'gcc' failed with exit status 1", one of the parameters is not assigned a value
# NOTE: if code yields error: "IndexError: index 0 is out of bounds for axis 0 with size 0", the location is likely wrong
# NOTE: if code yields warning: "ufunc 'isfinite' not supported for the input types" or "ValueError: cannot convert float NaN to integer",
#       population densities are very high, but this is not a problem for the density-independent model
#       if densities are not plotted, check the output .csv file and remove any rows with 0's after largest density
# NOTE: Ignoring "UserWarning: The target time is smaller than the current time" because the sampling step is small (which is not a problem
#       for jitcdde) and there is no backwards integration in time


# IMPORT PACKAGES
# Install and import jitcdde from https://github.com/neurophysik/jitcdde
import sys
import subprocess
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'jitcdde']) # check to see if jitcdde has been installed and if not, install it
from jitcdde import jitcdde, y, t

# Import other necessary packages
from numpy import arange, hstack, vstack, savetxt
from symengine import exp, pi, cos
from matplotlib.pyplot import subplots, xlabel, ylabel, xlim, ylim, yscale
from pandas import read_csv
from jitcxde_common import conditional
from os import chdir, getcwd

# Set working directory
chdir("..")


# USER: Enter species name and location or set all_pops to True
species = "Clavigralla shadabi"
location = "Benin"
all_pops = False

# USER: Run model for recent period (True) or future period (false)?
recent = True
if recent == True:
    period = "Recent"
else:
    period = "Future"

# USER: Save population dynamics data?
save = False

# USER: Include competition (i.e., density-dependent population growth)?
competition = True

# USER: Is model fit to census data?
census = False


# INPUT TEMPERATURE RESPONSE PARAMETERS AND HABITAT TEMPERATURE PARAMETERS
TR_params = read_csv(getcwd() + "/Model parameters/Temperature response parameters.csv")
T_params = read_csv(getcwd() + "/Model parameters/Habitat temperature parameters.csv")


# DEFINE MODEL PARAMETERS
# Time parameters
yr = 365 # days in year
if competition == True:
    start_yr = 0 # how many years before starting population dynamics model
    num_yrs = 75 # how long to run model
else: # NOTE: density-independent model is run between years 73 and 75 because otherwise densities become too large due to unbounded population growth
    start_yr = 73
    num_yrs = 2
if census == True:
    num_yrs = 10
    if species == "Clavigralla shadabi":
        start_yr = 1982
    elif species == "Apolygus lucorum":
        start_yr = 1996
        active_date = 137 # NOTE: for Apolygus lucorum during the census, insects only become active once the fields were planted on in mid-May (day = 137)
    
# Initial abundances
initJ = 100.
initA = 10.


# SELECT POPULATION
if all_pops == True:
    sp = 0
else:
    sp = TR_params[TR_params["Population"] == species + " " + location].index[0]


# START LOOP
while(True):
    
    # ASSIGN PARAMETERS
    sp_params = TR_params.iloc[sp]
    t_params = T_params.iloc[sp]
    
    # Set minimum tolerance to avoid integration errors associated with development rate approaching zero
    tol_min = 1e-8
    if recent == False and competition == True and (sp == 14 or sp == 19): # for Myzus persicae in Canada Chatham and US Columbia, minimum tolerance must be slightly higher
        tol_min = 1e-7
    elif recent == False and competition == True and sp == 17: # for Aulacorthum solani US Ithaca, minimum tolerance must be even higher
        tol_min = 1e-4
    
    # Temperature parameters (see "Habitat temperature parameters.csv")
    if recent == True:
        meanT = t_params["meanT.r"]
        amplT = t_params["amplT.r"] 
        shiftT = t_params["shiftT.r"]
        delta_mean = t_params["delta_mean.r"]
        delta_ampl = t_params["delta_ampl.r"]
    else:
        meanT = t_params["meanT.f"]
        amplT = t_params["amplT.f"] 
        shiftT = t_params["shiftT.f"]
        delta_mean = t_params["delta_mean.f"]
        delta_ampl = t_params["delta_ampl.f"]
    
    # Temperature responses (see "Temperature response parameters.csv")
    # per capita birth rate
    bTopt = sp_params["bTopt"]
    Toptb = sp_params["Toptb"]
    sb = sp_params["sb"]
    # development rate
    gMax = sp_params["gMax"]
    gTR = sp_params["gTR"]
    TR = sp_params["TR"]
    Ag = sp_params["Ag"]
    AL = sp_params["AL"]
    TL = sp_params["TL"]
    Tmin = sp_params["Tmin"]
    Toptg = sp_params["Toptg"]
    Tmaxg = sp_params["Tmaxg"]
    # stage-specific mortality rates
    dJTR = sp_params["dJTR"]
    AdJ = sp_params["AdJ"]
    dATR = sp_params["dATR"]
    AdA = sp_params["AdA"]
    # competition
    if competition == True:
        aMax = sp_params["aMax"]
    else:
        aMax = 0 # NOTE: can cause an "Unsuccessful Integration" error
    Topta = Toptb
    sa = sb
    
    
    # FUNCTIONS (where x is time because t is reserved for the state variable time in the DDE model)
    # Habitat temperature function (Eq. 5)
    def T(x):
            return conditional(x, 0, meanT, # during DDE "pre-history", habitat temperature is constant at its mean
                               (meanT + delta_mean*(x+start_yr*yr)) - (amplT + delta_ampl*(x+start_yr*yr)) * cos(2*pi*(x + shiftT)/yr))

    # Temperature response functions
    # per capita birth rate (Eq. 2a; see Eq. S3a for overwintering)
    if census == True and species == "Apolygus lucorum":
        def b(x):
            return conditional(T(x), T(active_date), 0, bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2))) # for the model assocaited with the Apolygus lucorum census, if habitat temperature is less than the temperature in mid-may when the field were planted, then the insect is overwintering
    else:
        def b(x):
            return conditional(T(x), Tmin, 0, bTopt * exp(-(T(x)-Toptb)**2/(2*sb**2))) # if habitat temperature < developmental min (Tmin), then use Eq. S3a; otherwise, use Eq. 2a
    
    # stage-specific per capita mortality rates (Eq. 2b; see Eq. S3b for overwintering)
    def dJ(x):
        return dJTR * exp(AdJ * (1/TR - 1/T(x)))
    def dA(x):
        return conditional(T(x), Tmin, 0.1 + dATR * exp(AdA * (1/TR - 1/T(x))), dATR * exp(AdA * (1/TR - 1/T(x)))) # if habitat temperature < developmental min (Tmin), then use Eq. S3b; otherwise, use Eq. 2b
   
    # development rate (Eq. 2c)
    def g(x):
        return conditional(gTR * T(x)/TR * exp(Ag * (1/TR - 1/T(x))) / (1 + exp(AL*(1/TL-1/T(x)))), tol_min, tol_min, # If g(T) < minimum tolerance, then g(T) = minimum tolerance (lines 96-100)
                           gTR * T(x)/TR * exp(Ag * (1/TR - 1/T(x))) / (1 + exp(AL*(1/TL-1/T(x))))) # If habitat temperature < developmental optima (Toptg), then use monotonic g(T)
                           #conditional(T(x), Toptg, gTR * T(x)/TR * exp(Ag * (1/TR - 1/T(x))) / (1 + exp(AL*(1/TL-1/T(x)))), # If habitat temperature < developmental optima (Toptg), then use monotonic g(T)
                                       #conditional(T(x), Tmaxg, gMax, tol_min)))  # If habitat temperature < developmental maximum (Tmaxg), then g(T) = gMax; otherwise, g(T) = minimum tolerance (lines 96-100)
        
    # density-dependence due to competition
    def a(x):
        return aMax * exp(-(T(x)-Topta)**2/(2*sa**2))
    
    
    # DDE MODEL
    # Define state variables
    J,A,S,τ = [y(i) for i in range(4)]
    
    # DDE model
    f = {
        J: b(t)*A*exp(-a(t)*A) - b(t-τ)*y(1,t-τ)*exp(-a(t-τ)*y(1,t-τ))*g(t)/g(t-τ)*S - dJ(t)*J, # juvenile density (Eq. 3a incorporating Eq. 4a); NOTE: python names this variable "y(0)"
        
        A: b(t-τ)*y(1,t-τ)*exp(-a(t-τ)*y(1,t-τ))*g(t)/g(t-τ)*S - dA(t)*A, # Adult density (Eq. 3b incorporating Eq. 4a); NOTE: python names this variable "y(1)"
        
        S: S*(g(t)/g(t-τ)*dJ(t-τ) - dJ(t)), # Juvenile through-stage survival (Eq. 4c); NOTE: python names this variable "y(2)"
        
        τ: 1 - g(t)/g(t-τ) # Developmental time-delay (Eq. 4b); NOTE: python names this variable "y(3)"
        }
    
    
    # RUN DDE SOLVER
    # Time and initial conditions
    times = arange(0, num_yrs*yr, 1)
    init = [ initJ, initA, exp(-dJ(-1e-3)/g(-1e-3)), 1./g(-1e-3)] # NOTE: using -1e-3 because 0 causes problems with model initialization
    
    # Run DDE solver
    DDE = jitcdde(f, max_delay=1e5, verbose=False)
    DDE.constant_past(init)
    DDE.compile_C(simplify=False, do_cse=False, verbose=True)
    DDE.adjust_diff()
    
    
    # SAVE DATA
    # array containing time and state variables
    data = vstack([ hstack([time, DDE.integrate(time)]) for time in times ])
    data[data < 1e-5] = 0 # set values below 1e-5 to 0
    
    # save data to csv
    if save == True:
        if census == False:
            if competition == True:
                filename = 'Time series data new/' + period + ' time series ' + sp_params["Population"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
            else:
                filename = 'Time series data DI new/' + period + ' time series ' + sp_params["Population"] + '.csv'
                savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
        else:
            filename = 'Time series data/Census time series ' + sp_params["Population"] + '.csv'
            savetxt(filename, data, fmt='%s', delimiter=",", header="Time,J,A,S,tau", comments='')
    
    
    # PLOT
    fig,ax = subplots()
    ax.plot(data[:,0], data[:,1], label='J') # plot juvenile density
    ax.plot(data[:,0], data[:,2], label='A') # plot adult density
    #ax.plot(data[:,0], vstack([T(i) - 273 for i in arange(0,num_yrs*yr,1) ]), label='T') # plot daily temperature
    #ax.plot(data[:,0], data[:,3], label='S') # plot juvenile survival
    #ax.plot(data[:,0], data[:,4], label='τ') # plot development time
    ax.legend(loc='best')
    xlabel("time (days)")
    ylabel("population density")
    yscale("linear")
    xlim((num_yrs-2)*yr,(num_yrs-0)*yr) # plot only last2 years of data
    ylim(0,100)
    
    
    # EXIT LOOP WHEN MODEL HAS BEEN RUN FOR SPECIFIED POPULATION OR ALL POPULATIONS
    if all_pops == True:
        sp = sp + 1
    else:
        sp = TR_params.shape[0] # if all_pops == False, then set "sp" to a value that will cause the loop to break
    if sp > TR_params.shape[0] - 1:
        break