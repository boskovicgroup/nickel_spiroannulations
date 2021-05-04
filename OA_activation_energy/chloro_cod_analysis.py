#!/usr/bin/env python
# coding: utf-8

# ## Activation energy from NMR spectra
# 0. Imports

# In[1]:


from kinetics.initial_rates import compInitialRate, compLeastSquaresRateLaw
from kinetics.json_definitions import readNMRIntegrationData

import os
import numpy
from matplotlib import pyplot as plt


# 1. Read the data from each NMR kinetic experiment JSON file.

# In[2]:


path_to_kinetic_data = 'chloro_substrate/cod/kinetic_json'

experiments = {}

# Iterate over the kinetic experiment directories.
for file_name in [f for f in os.listdir(path_to_kinetic_data) if not f.startswith('.')]:
    path_to_file = os.path.join(path_to_kinetic_data, file_name)
    
    # Read kinetic experiment data.
    kinetic_experiment = readNMRIntegrationData(path_to_file)
    
    # Record time points, integrals and temperature in the experiments dictionary
    experiment_id = kinetic_experiment["experiment_id"]
    experiments[experiment_id] = {
        "time_points": kinetic_experiment["time_points"],
        "integrals": kinetic_experiment["integrals"],
        "temperature": kinetic_experiment["temperature"]
    }


# 2. Compute the initial rates by estimating the initial rate of change in the integral value of the monitored component.

# In[3]:


experiments.keys()


# In[4]:


# Add reaction information associated to each component ordered according to its position in the component id list.
component_ids = ['component_1']
component_is_reactant = [False]
stoichiometric_coefficients = [1.0]

# Set parameters used to compute the initial rate of change in concentration of each component in each experiment.
experiments['cod_pos50']['initial_rate_parameters'] = {
    'component_1': {
        'time_range': slice(0,len(experiments['cod_pos50']['time_points'])),
        'reaction_percentage': 10
    },
}

experiments['cod_pos30']['initial_rate_parameters'] = {
    'component_1': {
        'time_range': slice(7,len(experiments['cod_pos30']['time_points'])),
        'reaction_percentage': 100
    },
}

experiments['cod_pos10']['initial_rate_parameters'] = {
    'component_1': {
        'time_range': slice(8,len(experiments['cod_pos10']['time_points'])),
        'reaction_percentage': 30
    },
}
experiments['cod_pos40']['initial_rate_parameters'] = {
    'component_1': {
        'time_range': slice(1,len(experiments['cod_pos40']['time_points'])),
        'reaction_percentage': 50
    },
}    
experiments['cod_pos20']['initial_rate_parameters'] = {
    'component_1': {
        'time_range': slice(3,len(experiments['cod_pos20']['time_points'])),
        'reaction_percentage': 50
    },
}

experiment_ids = list(experiments.keys())
for experiment_id in experiment_ids:
    initial_rates = []
    for component_index, component_id in enumerate(component_ids):
        experiment = experiments[experiment_id]
        time_range = experiment['initial_rate_parameters'][component_id]['time_range']
        reaction_percentage = experiment['initial_rate_parameters'][component_id]['reaction_percentage']
        time_points_subset = experiment['time_points'][time_range]
        integrals_subset = experiment['integrals'][time_range]
        stoichiometric_coefficient = stoichiometric_coefficients[component_index]
        t = time_points_subset[:]
        time_points_subset = list(numpy.array(time_points_subset)-time_points_subset[0])

        initial_rate = compInitialRate(stoichiometric_coefficient, time_points_subset, integrals_subset, reaction_percentage)
        initial_rates.append(initial_rate)
        
        # Generate initial rate curve.
        is_reactant = component_is_reactant[component_index]
        initial_rate_curve = integrals_subset[0]  + numpy.multiply(
            (-1 if is_reactant else 1) * stoichiometric_coefficient * initial_rate,
            time_points_subset)

        # Plot kinetic data and the initial rate curve.
        plt.title(f"{experiment_id}") 
        plt.xlabel("t/s") 
        plt.ylabel(f"Integral") 
        plt.plot(t, integrals_subset, "-bo", linestyle='none', label = f"Signal at t")
        plt.plot(t, initial_rate_curve,  "-r", label = "Initial rate")
        plt.legend(loc = "upper right")
        plt.ylim(min(integrals_subset), max(integrals_subset))
        plt.xlim(0, t[-1])
        plt.show()
        print("initial_rate = ", initial_rate, "signal/s")
    experiments[experiment_id]['initial_rates'] = initial_rates


# 3. Compute the activation energy using linear regression.

# In[5]:


# R = 8.31446261815324 # J/(mol*K)
R = 1.98720425864083e-3 # kcal/(mol*K)
initial_rates = [
    experiments[experiment_id]['initial_rates'][0]
    for experiment_id in experiment_ids
]
temperatures = [
    experiments[experiment_id]['temperature']
    for experiment_id in experiment_ids
]
A = numpy.vstack((numpy.ones(len(temperatures)), numpy.reciprocal(numpy.multiply(-R,temperatures)).T)).T
b = numpy.log(initial_rates)
x = numpy.linalg.lstsq(A, b, rcond=None)[0]
Ea = x[1]
fit = [x[0]-Ea/(R*temperature) for temperature in temperatures]
print(f"Ea: {x[1]} kcal/mol")
plt.xlabel("Temperature K") 
plt.ylabel("Log (Initial rate)") 
plt.plot(numpy.reciprocal(temperatures), numpy.log(initial_rates), "-bo", linestyle='none')
plt.plot(numpy.reciprocal(temperatures), fit, "-r")
plt.show()


# In[6]:


from matplotlib import rcParams
rcParams['font.sans-serif'] = "arial"
rcParams.update({'font.size': 22})
rcParams['legend.loc'] = 'upper right'
def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = numpy.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


# In[7]:


fig = plt.figure(figsize=(10, 10))
plt.plot(experiments['cod_pos10']['time_points'],experiments['cod_pos10']['integrals'], label = '10 \u00B0C', color = '#377eb8')
plt.plot(experiments['cod_pos20']['time_points'],experiments['cod_pos20']['integrals'], label = '20 \u00B0C', color = '#ff7f00')
plt.plot(experiments['cod_pos30']['time_points'],experiments['cod_pos30']['integrals'], label = '30 \u00B0C', color = '#4daf4a')
plt.plot(experiments['cod_pos40']['time_points'],experiments['cod_pos40']['integrals'], label = '40 \u00B0C', color = '#e41a1c')
plt.plot(experiments['cod_pos50']['time_points'],experiments['cod_pos50']['integrals'], label = '50 \u00B0C', color = '#984ea3')
plt.legend()
plt.xlabel('time [seconds]')
plt.ylabel('integral area')
plt.title('Formation of free COD; chloro substrate')
abline(3,0)
abline(10.42,0)
abline(17.27,0)
abline(88.10,0)
abline(305.70,0)
plt.xlim(0,5000)
plt.ylim(0,160000)
plt.savefig('/Users/zarko/Documents/zarkolab/projects/nickel_spiroindane/COD_formation_kinetics.svg')


# In[10]:


fig = plt.figure(figsize=(7, 6))
plt.plot(numpy.reciprocal(temperatures), numpy.log(initial_rates), marker = 'o', linestyle='none', color = '#377eb8')
plt.xlabel("1/temperature [1/K]") 
plt.ylabel("Log (initial rate)") 
plt.title('Arrhenius analysis; \n chloro COD formation')
plt.text(0.0033, 4, r'$E_a$ = 22.18 kcal/mol')
abline(-22.18203279/1.98720425864083e-3,40.11692398)
plt.savefig('/Users/zarko/Documents/zarkolab/projects/nickel_spiroindane/arrhenius_cl_cod.svg')


# In[ ]:




