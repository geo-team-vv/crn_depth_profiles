import os, sys
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats

# -----------------------------------------------------------------------------
# Function to list all the parameters of the CRN accumulation function into one container. 
# The output is a dictionary of values. 
def define_default_parameters():

    # Bulk density of overburden (2.1 per default)
    rho = 2.1

    # Production rate and production of 10Be (at/g/yr):
    # Prod rate follows Laloy et al. (2017), Sea Level High Latitude = 4.25, scaled after Stone (2000).
    prod_rate = 4.51

    # Production rate and production of 26Al:
    surf_prod_ratio = 6.75
    deep_prod_ratio = 7.4

    # Cosmogenic radionuclide constants and parameters.
    # atn=Attenuation length for neutrons (g/cm2), atnm=Attenuation length for negative muons, atfm=Attenuation length for fast muons
    atn = 152
    atnm = 1500
    atfm = 4320

    # Relative contribution to 10Be production:
    # pn, pnm, pfm= relative contribution of neutrons, negative muons and fast muons to 10Be production (unitless)
    pn = 0.9886
    pnm = 0.0027
    pfm = 0.0087

    # Relative contribution to 26Al production:
    # pn_Al, pnm_Al, pfm_Al = relative contribution of neutrons, negative muons and fast muons, respectively, to 26Al production (unitless)
    pn_Al = 0.9732
    pnm_Al = 0.0246
    pfm_Al = 0.0022

    # Half-life and decays for 10Be (yr) following Chmeleff (2010)
    half_life_Be10 = 1387000.0
    # Life and decays for 26Al (yr) following Nishiizumi (2004)
    half_life_Al = 705000.0

    # Time step (time resolution for aggradation phases in yr). 
    dt = 1000.0

    # Depth range, in cm, goes down from 0 cm at the surface to 4000 cm per default. 
    depth_values = range(4000)
    
    # initial_thickness is the thickness of the sediment column at t0 (cm). The erosion value at any time, and the difference between the total erosion and the total aggradation, must be lower than the base column thickness. Adapt the default value (1500 cm) if needed, otherwise the base column will become negative and make the model crash.
    initial_thickness = 1500

    # Dictionary of parameters
    parameters = {
        "rho": rho,
        "prod_rate": prod_rate,
        "surf_prod_ratio": surf_prod_ratio,
        "deep_prod_ratio": deep_prod_ratio,
        "atn": atn,
        "atnm": atnm,
        "atfm": atfm,
        "pn": pn,
        "pnm": pnm,
        "pfm": pfm,
        "pn_Al": pn_Al,
        "pnm_Al": pnm_Al,
        "pfm_Al": pfm_Al,
        "half_life_Be10": half_life_Be10,
        "half_life_Al": half_life_Al,
        "dt": dt,
        "depth_values": depth_values,
        "initial_thickness": initial_thickness,
    }
    
    return parameters

# -----------------------------------------------------------------------------
# Function to compute the nuclide production rate at a given depth. 
def prod_rate_depth(depth, prod_rate_surf, rho, attenuation):
     prod_depth = prod_rate_surf*math.exp((-depth*rho) / attenuation)

     return prod_depth

# -----------------------------------------------------------------------------
# Fonction to compute the nuclide concentration at a given depth
def conc_depth(N_acc, prod_rate_depth, decay, timestep, duration_remaining, erosion_rate, rho, attenuation):

    conc_depth = N_acc + (prod_rate_depth / (decay + (erosion_rate*rho)/attenuation)) * (1 - math.exp(-timestep*(decay + (erosion_rate*rho)/attenuation)))*math.exp(-decay*duration_remaining)

    return conc_depth

# -----------------------------------------------------------------------------
# Function to compute the Reduced chi-squared value
def compute_reduced_chisquare(obs_data_Be, obs_data_Al, fitted_concentrations, degrees_freedom = 10):

    # Merge 26Al and 10Be observed data into arrays
    conc_observed = np.concatenate((obs_data_Be["observed_concentration"].to_numpy(), obs_data_Al["observed_concentration"].to_numpy()))
    sd_observed = np.concatenate((obs_data_Be["observed_sd"].to_numpy(), obs_data_Al["observed_sd"].to_numpy()))

    # Concatenated fitted Al and Be concentrations values
    conc_fitted = np.concatenate((fitted_concentrations["fitted_conc_Be"], fitted_concentrations["fitted_conc_Al"]))

    if sd_observed[0] == None:
       chisq = np.sum((conc_observed - conc_fitted)**2)
    else:
       chisq = np.sum(((conc_observed - conc_fitted) / sd_observed)**2)

    nu = conc_observed.size - 1 - degrees_freedom

    chi_value = chisq / nu
    p_value = 1 - stats.chi2.cdf(chi_value*degrees_freedom, degrees_freedom)

    chi = {
        "chi_value": chi_value,
        "p_value": p_value,
    }

    return chi

# -----------------------------------------------------------------------------
# Function to get fitted values at observed depths
def get_fitted_concentrations(data_Be, data_Al, depth_profiles):

    # Containers for CRN concentration fitted values
    fitted_conc_Be = []
    fitted_conc_Al = []

    # Extract fitted CRN concentrations at sampled depths
    for i in range(len(data_Be["depth"])):
        fitted_conc_Be.append(depth_profiles["crn_conc_Be"][(data_Be["depth"][i])])

    for i in range(len(data_Al["depth"])):
        fitted_conc_Al.append(depth_profiles["crn_conc_Al"][data_Al["depth"][i]])

    # Return fitted CRN concentrations every centimeter for Al and Be
    fitted_concentrations = {
        "fitted_conc_Be": np.array(fitted_conc_Be),
        "fitted_conc_Al": np.array(fitted_conc_Al),
    }

    return fitted_concentrations

# -----------------------------------------------------------------------------
# CRN accumulation over depositional and post depositional history
# total_time, geomorpho_history and N_inh are independant lists (type values between brackets, as for classical lists).
# total_time is a list of time periods in years. exposure_age is the sum of all the time periods specified in the list t_Fin.
# geomorpho_history is a list of denudation or aggradation values in cm. By convention, erosion is positive and aggradation negative. A value of 0 means landscape stability. In case of aggradation period, the amount of centimeters the user wishes to aggrade must be preceded by a minus sign. total_time and Erosion_Burial lists MUST have the same length.
# N_inh is a list of inherited 10Be concentrations. The first term of the list is the inherited concentration of the material already deposited at the beginning of the simulation. The vector has the same length as the Total_time vector. For each parametrized erosion or stability period, just insert a N_inh term =0. 
# The other terms in N_inh are the inherited concentrations of the material accumulated during each aggradation phas. N_inh must thus be written with zero values when positions correspond to an erosion or a stability phase, and with an inherited concentration value when the position corresponds to an aggradation phase.
# exposure age remaining first equals exposure age at t0 and will decrease as we move from one period to another. 

def compute_crn_depth_profile(total_time, geomorpho_history, N_inh, parameters):
    
    # -----------------------------------------------------------------------------
    # Compute secondary variables
    # prod_n, prod_nm and prod_fm = fraction of the production rate respectively provided by spallation, negative muons and fast muons (at/g/yr)
    prod_n_Be = parameters["prod_rate"] * parameters["pn"]
    prod_nm_Be = parameters["prod_rate"] * parameters["pnm"]
    prod_fm_Be = parameters["prod_rate"] * parameters["pfm"]

    # L_Be = decay constant (lambda) for 10Be
    L_Be = math.log(2) / parameters["half_life_Be10"]

    # Production rate and production of 26Al:
    # surf_prod_ratio = ratio of production between 26Al and 10Be (surface), deep_prod_ratio = mean between the surface production ratio and the production between 26Al and 10Be at great depth (variable between studies, 8.0 for Margreth et al. (2016) and Knudsen et al. (2019), 8.4 for Akcar et al., 2017).
    # exposure_Al = Production rate of 10Be times surface production ratio. prod_n_Al, prod_nm_Al, prod_fm_Al= fraction of the 26Al production rate respectively provided by spallation, negative muons and fast muons (at/g/yr)

    exposure_Al = parameters["prod_rate"] * parameters["surf_prod_ratio"]
    prod_n_Al = exposure_Al * parameters["pn_Al"]
    prod_nm_Al = exposure_Al * parameters["pnm_Al"]
    prod_fm_Al = exposure_Al * parameters["pfm_Al"]

    ## L_Al= decay constant (lambda) for 26Al
    L_Al = math.log(2) / parameters["half_life_Al"]

    # -----------------------------------------------------------------------------
    # Generate depth profiles of production rates. 
    # prod_depth_spal_Be, prod_depth_neg_muon_Be, prod_depth_fast_muon_Be = 10Be production distributed over depth from spallation, negative muons and fast muons, respectively.
    # prod_depth_spal_Al, prod_depth_neg_muon_Al, prod_depth_fast_muon_Al = 26Al production distributed over depth from spallation, negative muons and fast muons, respectively.
    prod_depth_spal_Be = np.vectorize(prod_rate_depth)(parameters["depth_values"], prod_n_Be, parameters["rho"], parameters["atn"])
    prod_depth_neg_muon_Be = np.vectorize(prod_rate_depth)(parameters["depth_values"], prod_nm_Be, parameters["rho"], parameters["atnm"])
    prod_depth_fast_muon_Be = np.vectorize(prod_rate_depth)(parameters["depth_values"], prod_fm_Be, parameters["rho"], parameters["atfm"])
    prod_depth_spal_Al = np.vectorize(prod_rate_depth)(parameters["depth_values"], prod_n_Al, parameters["rho"], parameters["atn"])
    prod_depth_neg_muon_Al = np.vectorize(prod_rate_depth)(parameters["depth_values"], prod_nm_Al, parameters["rho"], parameters["atnm"])
    prod_depth_fast_muon_Al = np.vectorize(prod_rate_depth)(parameters["depth_values"], prod_fm_Al, parameters["rho"], parameters["atfm"])

    exposure_age = sum(total_time)
    exposure_age_remaining = exposure_age

    # 10Be Column creation:
    # crn_conc_Be are empty arrays building a 1500 cm depth space where to store the values of accumulated 10Be from spallation, negative muons and fast muons
    crn_conc_spal_Be = np.zeros(parameters["initial_thickness"])
    crn_conc_neg_muon_Be = np.zeros(parameters["initial_thickness"])
    crn_conc_fast_muon_Be = np.zeros(parameters["initial_thickness"])
    crn_conc_Be = np.zeros(parameters["initial_thickness"])

    # 26Al Column creation:
    # crn_conc_Be are empty arrays building a 1500 cm depth space where to store the values of accumulated 26Al from spallation, negative muons and fast muons
    crn_conc_spal_Al = np.zeros(parameters["initial_thickness"])
    crn_conc_neg_muon_Al = np.zeros(parameters["initial_thickness"])
    crn_conc_fast_muon_Al = np.zeros(parameters["initial_thickness"])
    crn_conc_Al = np.zeros(parameters["initial_thickness"])

    # Steps and indices:
    # k = 0 is the onset of index to move from one period and associated geomorpho_history value to another.
    k = 0

    # Be and Al Inheritance and decay in the base column (at t = 0)
    # N_inh_fin_Be is the present day 10Be inheritance. It is the inheritance at t = 0 (first term of the list of parameters N_inh) that has decayed over the exposure age (sum of total_time parts). This is stored in crn_conc_Be that collects the final concentrations from all the production paths and inheritance. Regarding 26Al, inheritance in the column present at t = 0 assumes that all the deposits got a surface inherited ratio (6.75 times the amount of inherited 10Be). You may replace the parameter parameters["surf_prod_ratio"] by a value of your own if you estimate this ratio different.
    N_inh_fin_Be = N_inh[0]*math.exp(-(L_Be*exposure_age))
    N_inh_fin_Al = N_inh[0]*parameters["surf_prod_ratio"]*math.exp(-(L_Al*exposure_age))
    
    # for i in range(len(crn_conc_Be)):
    crn_conc_Be = crn_conc_Be + N_inh_fin_Be
    
    # for i in range(len(crn_conc_Al)):
    crn_conc_Al = crn_conc_Al + N_inh_fin_Al

    # CRN Accumulation during the different time period parametrized in list total_time
    for k in range(len(total_time)):

        # Stability period
        if geomorpho_history[k] == 0.0:

            exposure_age_remaining -= total_time[k]
            
            if exposure_age_remaining < 0.0:
                exposure_age_remaining = 0.0
            
            # Get current column depth
            current_column_depth = len(crn_conc_spal_Be)

            # For spallation, compute CRN concentration over the depth column.
            crn_conc_spal_Be = conc_depth(N_acc = crn_conc_spal_Be, prod_rate_depth = prod_depth_spal_Be[:current_column_depth], decay = L_Be, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atn"])
            crn_conc_spal_Al = conc_depth(N_acc = crn_conc_spal_Al, prod_rate_depth = prod_depth_spal_Al[:current_column_depth], decay = L_Al, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atn"])
            # For negative muons, compute CRN concentration over the depth column.
            crn_conc_neg_muon_Be = conc_depth(N_acc = crn_conc_neg_muon_Be, prod_rate_depth = prod_depth_neg_muon_Be[:current_column_depth], decay = L_Be, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atnm"])
            crn_conc_neg_muon_Al = conc_depth(N_acc = crn_conc_neg_muon_Al, prod_rate_depth = prod_depth_neg_muon_Al[:current_column_depth], decay = L_Al, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atnm"])
            # For fast muons, compute CRN concentration over the depth column.
            crn_conc_fast_muon_Be = conc_depth(N_acc = crn_conc_fast_muon_Be, prod_rate_depth = prod_depth_fast_muon_Be[:current_column_depth], decay = L_Be, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atfm"])
            crn_conc_fast_muon_Al = conc_depth(N_acc = crn_conc_fast_muon_Al, prod_rate_depth = prod_depth_fast_muon_Al[:current_column_depth], decay = L_Al, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atfm"])

        # Erosion period
        # erosion_rate is the erosion rate. eroded_thickness_cm is the amount of top centimeters removed due to erosion, eroded_thickness_cm = int(eroded_thickness_cm) rounds up the number of centimeter eroded to the closest integer value.
        elif geomorpho_history[k] > 0.0:

            exposure_age_remaining -= total_time[k]
            erosion_rate = geomorpho_history[k] / total_time[k]
            
            if exposure_age_remaining < 0.0:
                exposure_age_remaining = 0.0

            eroded_thickness_cm = geomorpho_history[k]
            eroded_thickness_cm = int(eroded_thickness_cm)

            crn_conc_spal_Be = np.delete(crn_conc_spal_Be, range(eroded_thickness_cm))
            crn_conc_neg_muon_Be = np.delete(crn_conc_neg_muon_Be, range(eroded_thickness_cm))
            crn_conc_fast_muon_Be = np.delete(crn_conc_fast_muon_Be, range(eroded_thickness_cm))
            crn_conc_Be = np.delete(crn_conc_Be, range(eroded_thickness_cm))
            crn_conc_spal_Al = np.delete(crn_conc_spal_Al, range(eroded_thickness_cm))
            crn_conc_neg_muon_Al = np.delete(crn_conc_neg_muon_Al, range(eroded_thickness_cm))
            crn_conc_fast_muon_Al = np.delete(crn_conc_fast_muon_Al, range(eroded_thickness_cm))
            crn_conc_Al = np.delete(crn_conc_Al, range(eroded_thickness_cm))

            # Get current column depth
            current_column_depth = len(crn_conc_spal_Be)

            # For spallation, compute CRN concentration over the depth column.
            crn_conc_spal_Be = conc_depth(N_acc = crn_conc_spal_Be, prod_rate_depth = prod_depth_spal_Be[:current_column_depth], decay = L_Be, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = erosion_rate, rho = parameters["rho"], attenuation = parameters["atn"])
            crn_conc_spal_Al = conc_depth(N_acc = crn_conc_spal_Al, prod_rate_depth = prod_depth_spal_Al[:current_column_depth], decay = L_Al, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = erosion_rate, rho = parameters["rho"], attenuation = parameters["atn"])
            # For negative muons, compute CRN concentration over the depth column.
            crn_conc_neg_muon_Be = conc_depth(N_acc = crn_conc_neg_muon_Be, prod_rate_depth = prod_depth_neg_muon_Be[:current_column_depth], decay = L_Be, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = erosion_rate, rho = parameters["rho"], attenuation = parameters["atnm"])
            crn_conc_neg_muon_Al = conc_depth(N_acc = crn_conc_neg_muon_Al, prod_rate_depth = prod_depth_neg_muon_Al[:current_column_depth], decay = L_Al, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = erosion_rate, rho = parameters["rho"], attenuation = parameters["atnm"])
            # For fast muons, compute CRN concentration over the depth column.
            crn_conc_fast_muon_Be = conc_depth(N_acc = crn_conc_fast_muon_Be, prod_rate_depth = prod_depth_fast_muon_Be[:current_column_depth], decay = L_Be, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = erosion_rate, rho = parameters["rho"], attenuation = parameters["atfm"])
            crn_conc_fast_muon_Al = conc_depth(N_acc = crn_conc_fast_muon_Al, prod_rate_depth = prod_depth_fast_muon_Al[:current_column_depth], decay = L_Al, timestep = total_time[k], duration_remaining = exposure_age_remaining, erosion_rate = erosion_rate, rho = parameters["rho"], attenuation = parameters["atfm"])

        # Aggradation period
        # "exposure_age_remaining -= dt" means the removal of a first time step from exposure age to enter the accumulation loop over that time step. duration_period_k = 0 is the onset of the aggradation time period. aggrad_accumulated = 0 will store aggraded material (useful if centimetric values are not reached per time step, e.g. for very slow aggradation scenarios). aggradation_rate = geomorpho_history[k]/total_time[k] calculates the aggradation rate in cm/yr. aggraded_thickness_per_timestep_cm = aggradation_rate \_dt*(-1) calculates how much material is delivered as aggradation layers in cm per timestep. Aggradation is negative by convention in the model, so we apply \*(-1). duration_period_k+=parameters["dt"] means that the duration period is increased by one time step at every aggradation loop until it equals total_time[k]
        elif geomorpho_history[k] < 0.0:
            exposure_age_remaining -= parameters["dt"]
            duration_period_k = 0
            aggrad_accumulated = 0
            aggradation_rate = geomorpho_history[k] / total_time[k]
            aggraded_thickness_per_timestep_cm = aggradation_rate*parameters["dt"]*(-1)

            while duration_period_k < total_time[k]:
                duration_period_k += parameters["dt"]
                aggrad_accumulated += aggraded_thickness_per_timestep_cm

                # Everytime aggrad_accumulated is at least 1cm, this cm will be added on top of the column, with its related inheritance.this layer will then start to accumulate CRN via spallation, negative muons and fest muons until the end of the exposure time remaining
                if aggrad_accumulated >= 1.0:
                    while aggrad_accumulated >= 1.0:
                        aggrad_accumulated -= 1.0
                        crn_conc_spal_Be = np.insert(crn_conc_spal_Be, 0, 0.0)
                        crn_conc_neg_muon_Be = np.insert(crn_conc_neg_muon_Be, 0, 0.0)
                        crn_conc_fast_muon_Be = np.insert(crn_conc_fast_muon_Be, 0, 0.0)
                        crn_conc_spal_Al = np.insert(crn_conc_spal_Al, 0, 0.0)
                        crn_conc_neg_muon_Al = np.insert(crn_conc_neg_muon_Al, 0, 0.0)
                        crn_conc_fast_muon_Al = np.insert(crn_conc_fast_muon_Al, 0, 0.0)

                        N_inh_current = N_inh[k]*math.exp(-L_Be*exposure_age_remaining)
                        crn_conc_Be = np.insert(crn_conc_Be, 0, N_inh_current)

                        # N_inh_current_Al is multiplied by the surface production ratio, but can be multiplied by a deeper prodcution ratio (> 6.75) is the user wishes. In such case, simply replace the "parameters["surf_prod_ratio"]" in the next equation by the "deep_prod_ratio" (defined above in the code as being 7.4), or replace it directly by a value of your own.
                        N_inh_current_Al = N_inh[k]*parameters["surf_prod_ratio"]*math.exp(-L_Al*exposure_age_remaining)
                        crn_conc_Al = np.insert(crn_conc_Al, 0, N_inh_current_Al)

                # Get current column depth
                current_column_depth = len(crn_conc_spal_Be)

                # For spallation, compute CRN concentration over the depth column.
                crn_conc_spal_Be = conc_depth(N_acc = crn_conc_spal_Be, prod_rate_depth = prod_depth_spal_Be[:current_column_depth], decay = L_Be, timestep = parameters["dt"], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atn"])
                crn_conc_spal_Al = conc_depth(N_acc = crn_conc_spal_Al, prod_rate_depth = prod_depth_spal_Al[:current_column_depth], decay = L_Al, timestep = parameters["dt"], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atn"])
                # For negative muons, compute CRN concentration over the depth column.
                crn_conc_neg_muon_Be = conc_depth(N_acc = crn_conc_neg_muon_Be, prod_rate_depth = prod_depth_neg_muon_Be[:current_column_depth], decay = L_Be, timestep = parameters["dt"], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atnm"])
                crn_conc_neg_muon_Al = conc_depth(N_acc = crn_conc_neg_muon_Al, prod_rate_depth = prod_depth_neg_muon_Al[:current_column_depth], decay = L_Al, timestep = parameters["dt"], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atnm"])
                # For fast muons, compute CRN concentration over the depth column.
                crn_conc_fast_muon_Be = conc_depth(N_acc = crn_conc_fast_muon_Be, prod_rate_depth = prod_depth_fast_muon_Be[:current_column_depth], decay = L_Be, timestep = parameters["dt"], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atfm"])
                crn_conc_fast_muon_Al = conc_depth(N_acc = crn_conc_fast_muon_Al, prod_rate_depth = prod_depth_fast_muon_Al[:current_column_depth], decay = L_Al, timestep = parameters["dt"], duration_remaining = exposure_age_remaining, erosion_rate = 0, rho = parameters["rho"], attenuation = parameters["atfm"])

                exposure_age_remaining -= parameters["dt"]

                if exposure_age_remaining < 0.0:
                    exposure_age_remaining = 0.0

    # Sum up the different contributors to 26Al and 10Be production:
    crn_conc_Be = crn_conc_Be + crn_conc_spal_Be + crn_conc_neg_muon_Be + crn_conc_fast_muon_Be
    crn_conc_Al = crn_conc_Al + crn_conc_spal_Al + crn_conc_neg_muon_Al + crn_conc_fast_muon_Al

    depth_profiles = {
        "crn_conc_Be": crn_conc_Be,
        "crn_conc_Al": crn_conc_Al,
    }

    # Output: Depth profiles of CRN concentrations for 10Be and 26Al.
    return depth_profiles 

# -----------------------------------------------------------------------------
# Function that create arrays of Monte Carlo parameters from csv file
def define_monte_carlo_parameters(data_monte_carlo, n_simus):

    # Initialize containers for Monte-Carlo values for n_simus simulations
    total_times = np.empty((0, n_simus), int)
    geomorpho_histories = np.empty((0, n_simus), int)
    N_inh_values = np.empty((0, n_simus), int)

    # Loop within rows, i.e. steps of aggradation, erosion or stability periods and get min and max values
    for index, row in data_monte_carlo.iterrows():

        total_times = np.append(total_times, [np.rint(np.random.uniform(row["total_time_min"], row["total_time_max"], n_simus)).astype(int)], axis = 0)
        geomorpho_histories = np.append(geomorpho_histories, [np.rint(np.random.uniform(row["geomorpho_history_min"], row["geomorpho_history_max"], n_simus)).astype(int)], axis = 0)
        N_inh_values = np.append(N_inh_values, [np.rint(np.random.uniform(row["n_inh_min"], row["n_inh_max"], n_simus)).astype(int)], axis = 0)

    # Create dictionary containing Monte Carlo parameters
    monte_carlo_parameters = {
        "total_times": total_times,
        "geomorpho_histories": geomorpho_histories,
        "N_inh_values": N_inh_values
    }

    return monte_carlo_parameters

# -----------------------------------------------------------------------------
# Compute Monte Carlo simulations based on ranges of parameters
def compute_monte_carlo_simulations(params_monte_carlo, data_Be, data_Al, parameters):

    # List container for scenarios
    scenarios = []

    # Get number of simulations to compute
    n_simus = len(params_monte_carlo["total_times"][0])

    # Loop through number of simulations
    for simu in range(0, n_simus):

        # params as list 
        total_time = [values[simu] for values in params_monte_carlo["total_times"]]
        geomorpho_history = [values[simu] for values in params_monte_carlo["geomorpho_histories"]]
        N_inh = [values[simu] for values in params_monte_carlo["N_inh_values"]]

        # Compute modelized depth profile 
        depth_profiles = compute_crn_depth_profile(total_time, geomorpho_history , N_inh, parameters)

        # Get fitted concentrations at sampled depths
        fitted_concentrations = get_fitted_concentrations(data_Be, data_Al, depth_profiles)

        # Reduced chi-squared calculation
        chi = compute_reduced_chisquare(data_Be, data_Al, fitted_concentrations, degrees_freedom = 10)

        # Save scenario parameters as a possible solution if reduced chi-square is significant
        if chi["p_value"] >= 0.05:

            # Collect data into a dictionary
            current_scenario = {
                "total_time": total_time,
                "geomorpho_history": geomorpho_history,
                "N_inh": N_inh,
                "chi": chi
            }
            # Append current scenario to container
            scenarios.append(current_scenario)

            # Free memory from current objects
            del current_scenario

        # Free memory from current objects
        del depth_profiles, fitted_concentrations, chi

        # Print progress every 10% of simulations
        if simu % (n_simus/10) == 0:

            # Print every 10% of simulations
            print("## COMPUTING MONTE CARLO SIMULATIONS: ", simu/n_simus*100, "%")
    
    # Print number of scenarios
    print("## NUMBER OF SIGNIFICANT SCENARIOS (p-value >= 0.05): ", len(scenarios))

    # Return scenarios
    return scenarios

# -----------------------------------------------------------------------------
# Get variables from scenarios to be plotted for a given step
def get_variables_from_scenarios(scenarios, step):
    
    # Container for values extracted from scenarios to be plotted
    data_plot = {
        "total_time": [],
        "geomorpho_history": [],
        "chi_values": []
        }

    # Loop through stored scenarios
    for scenario in scenarios:

        # Collect values
        data_plot["total_time"].append(scenario["total_time"][step-1])
        data_plot["geomorpho_history"].append(scenario["geomorpho_history"][step-1])
        data_plot["chi_values"].append(scenario["chi"]["chi_value"])
    
    data_plot["step"] = step

    # Count number of scenarios for current step and print message
    if len(data_plot["total_time"]) == 0:
        print("## NO SIGNIFICANT SCENARIO FOR STEP HAVE BEEN MODELLED (data_plot is empty), i.e. you should increase the number of simulations (n_simus).")

    # Return data from scenarios to be plotted
    return data_plot

# -----------------------------------------------------------------------------
# Plot Gaussian Kernels 
def plot_gaussian_kernel(data_plot, data_monte_carlo, dir_root):

    # Check if data_plot contains values 
    if len(data_plot["total_time"]) == 0:
        print("## NO GRAPH CAN BE PLOTTED BECAUSE DATA_PLOT IS EMPTY.")
    
    # Plot graphs if data_plot is not empty
    else:

        # Create output directory if it does not exists
        dir_graphs = os.path.join(dir_root, "graphs")
        if not os.path.exists(dir_graphs):
            os.makedirs(dir_graphs)

        # Format the colorbar
        def fmt(x, pos): 
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

        step = data_plot["step"]

        # Get ranges of parameters for current step
        if data_monte_carlo.iloc[step - 1]["geomorpho_history_min"] < 0:
            xmin = data_monte_carlo.iloc[step - 1]["geomorpho_history_min"]
        else: 
            xmin = 0.0
        xmax = data_monte_carlo.iloc[step - 1]["geomorpho_history_max"]
        ymin = 0.0
        ymax = data_monte_carlo.iloc[step - 1]["total_time_max"]

        # Create meshgrid
        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

        # We will fit a gaussian kernel using the scipy’s gaussian_kde method:
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([np.array(data_plot["geomorpho_history"]), np.array(data_plot["total_time"])])
        kernel = stats.gaussian_kde(values, bw_method = 'scott')
        f = np.reshape(kernel(positions).T, xx.shape)

        # Plotting gauss kernel with filled contours
        plt.figure(figsize=(5, 4), facecolor = 'white')
        cfset = plt.contourf(xx, yy, f, cmap = 'jet') 
        plt.xlabel('U' + str(step) + ' erosion (cm)')
        plt.ylabel('U' + str(step) + ' hiatus (yr)')
        plt.axis([xmin, xmax, ymin, ymax])
        # 0.0024 = max(kernel)
        plt.colorbar(format = ticker.FuncFormatter(fmt)) 
        # plt.colorbar(ticks = [0.000, 0.0012, 0.0024], format = ticker.FuncFormatter(fmt)) 
        plt.savefig(os.path.join(dir_graphs, "graph-geom_hist-duration-step-" + str(step) +".png"))
