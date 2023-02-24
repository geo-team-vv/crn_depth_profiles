# Libraries to import
import timeit
from utilities.CRN_functions import *

# -----------------------------------------------------------------------------
# Read data files

# Get root directory
dir_root = os.path.dirname(sys.path[0])
# Read csv files containing observed data
data_Al = pd.read_csv(os.path.join(dir_root, "data", "obs_data_al.csv"), sep = ";", header = 0)
data_Be = pd.read_csv(os.path.join(dir_root, "data", "obs_data_be.csv"), sep = ";", header = 0)
# Read parameters values for Monte Carlo simulations
data_monte_carlo = pd.read_csv(os.path.join(dir_root, "data", "params_monte_carlo.csv"), sep = ";")

# -----------------------------------------------------------------------------
# Define parameters of CRN modelling
parameters = define_default_parameters()

# -----------------------------------------------------------------------------
# Simple Monte Carlo simulations

# Number of simulations
n_simus = 1000

# If true, the input parameters of Monte Carlo simulations, i.e. steps defined in the params_monte_carlo.csv
# file, will be modified according to Vandermaelen et al. (2022).
reproduce_scenario_geosite_vandermaelen = True

# Compute Monte Carlo parameters based on parameters file
params_mc = define_monte_carlo_parameters(data_monte_carlo, n_simus)

# -----------------------------------------------------------------------------
# Modify parameters that are functions of others, in order to constrain scenarios
# as they are presented in Vandermaelen et al. (2022), if parameter scenario_geosite_as is true.

if(reproduce_scenario_geosite_vandermaelen):

    # Compute total aggradation time without step 5, which is constrained by external data
    aggradation_duration_without_step_5 = params_mc["total_times"][0] + params_mc["total_times"][1] + params_mc["total_times"][2] + params_mc["total_times"][3] + params_mc["total_times"][5] + params_mc["total_times"][6] + params_mc["total_times"][7]
    params_mc["total_times"][4] = np.rint(np.random.uniform(500000.0, 1000000.0 - aggradation_duration_without_step_5, n_simus)).astype(int)

    # Define a normal distribution for inheritance of the first step and use these values for steps 2 and 4
    params_mc["N_inh_values"][0] = np.absolute(np.rint(np.random.normal(900000.0, 20000.0, n_simus)).astype(int))
    params_mc["N_inh_values"][1] = params_mc["N_inh_values"][0]
    params_mc["N_inh_values"][3] = params_mc["N_inh_values"][0]

    # Compute erosion values of steps 3, 6 and 8, which depends on previous aggradation values
    params_mc["geomorpho_histories"][2] = np.absolute(params_mc["geomorpho_histories"][1]) - 185
    params_mc["geomorpho_histories"][5] = np.absolute(params_mc["geomorpho_histories"][3]) - 275
    params_mc["geomorpho_histories"][7] = np.absolute(params_mc["geomorpho_histories"][6]) - 60

# -----------------------------------------------------------------------------
# Compute simple Monte Carlo simulations for n_simus simulations

# Track processing time
start_time = timeit.default_timer()
scenarios = compute_monte_carlo_simulations(params_mc, data_Be, data_Al, parameters)

# Print  processing time
print("TIME ELAPSED: " + str(timeit.default_timer() - start_time))

# -----------------------------------------------------------------------------
# Define step of history for which data will be extracted and plotted (Steps 1, 3 and 6 are relevant
# in the frame of Vandermaelen et al., 2022).
step = 6

# Get variables from scenarios to be plotted for a given step
data_plot = get_variables_from_scenarios(scenarios, step)

# -----------------------------------------------------------------------------
# Plot Gaussian Kernels 
plot_gaussian_kernel(data_plot, data_monte_carlo, dir_root)
