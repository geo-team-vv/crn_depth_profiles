# -----------------------------------------------------------------------------
# Libraries to import
from CRN_functions import *

# -----------------------------------------------------------------------------
# Read data files

# Get root directory
dir_root = os.path.dirname(sys.path[0])
# Read csv files containing observed data
data_Al = pd.read_csv(dir_root + "/data/obs_data_al.csv", sep = ";", header = 0)
data_Be = pd.read_csv(dir_root + "/data/obs_data_be.csv", sep = ";", header = 0)

# -----------------------------------------------------------------------------
# Define parameters of CRN modelling
parameters = define_parameters()

# -----------------------------------------------------------------------------
#  Implementation of one single scenario of erosion history

# Fill the total_time, geomorpho_history and N_inh list with your own scenario values
total_time = [10000.0, 10000.0, 10000.0, 10000.0, 10000.0]
geomorpho_history = [100, -50, 200, 0, -100]
N_inh = [10000, 0, 0, 0, 5000]

# -----------------------------------------------------------------------------
# Compute CRN accumulation along a depth profile
depth_profiles = compute_crn_depth_profile(total_time, geomorpho_history, N_inh, parameters) 
# Get fitted concentrations at sampled depths
fitted_concentrations = get_fitted_concentrations(data_Be, data_Al, depth_profiles)
print(fitted_concentrations)

# Reduced chi-squared calculation
chi = compute_reduced_chisquare(data_Be, data_Al, fitted_concentrations, degrees_freedom = 10)
print(chi)
