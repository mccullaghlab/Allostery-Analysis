
# ----------------------------------------
# USAGE AND REFERENCING
# ----------------------------------------
#   python create_hessian.py create_hessian_config_file_name

# ----------------------------------------
# CODE OUTLINE
# ----------------------------------------
#   Allosteric Paths in Proteins - Hessian Analysis Code
#       1) Load in user defined parameters
#       2) Load in necessary functions from module files
#       3) Hessian Analysis

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------
import sys
import os
import importlib
import numpy as np
from IO import create_hessian_config_parser,create_hessian_summary

# ----------------------------------------
# VARIABLE DECLARATION: 
# ----------------------------------------
config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------

def main():
        covariance = np.loadtxt(parameters['system_covariance_file'])
        average_structure = np.loadtxt(parameters['system_average_structure_file'])
        if parameters['initial_guess_hessian'] != None:
                guess_hessian = np.loadtxt(parameters['initial_guess_hessian'])
        else:
                guess_hessian = None

        # ----------------------------------------
        # 3) Hessian Analysis
        # ----------------------------------------
        if parameters['hessian_algorithm'].lower() == 'henm':
                henm_hessian = calculate_hessian(covariance, average_structure, parameters['output_directory'], guess = guess_hessian, max_iterations = int(parameters['henm_max_iterations']), alpha = float(parameters['henm_alpha']), kBT = 1.9872E-3*float(parameters['temperature']), threshold = float(parameters['henm_threshold']),plotting_boolean = parameters['plotting_boolean'],output_step=1,distance_cutoff = parameters['distance_cutoff'])
                np.savetxt(parameters['output_directory'] + parameters['output_hessian_file_name'],henm_hessian)
        #elif parameters['hessian_algorithm'].lower() == 'reach':
        #        reach_hessian = calculate_hessian()
        #        np.savetxt(parameters['output_directory'] + parameters['output_hessian_file_name'],reach_hessian)

        if parameters['summary_boolean']:
                create_hessian_summary(parameters['output_directory'] + 'create_hessian.summary',sys.argv,parameters)

# ----------------------------------------
# 1) LOAD IN USER DEFINED PARAMETERS
# ----------------------------------------
parameters = {}
create_hessian_config_parser(config_file,parameters)

# ----------------------------------------
# SETTING UP THE OUTPUT DIRECTORY
# ----------------------------------------
if parameters['output_directory'][-1] != os.sep:
        parameters['output_directory'] += os.sep

#if os.path.exists(parameters['output_directory']):
#        print 'The output directory, ', parameters['output_directory'], 'already exists. Please select a different directory name for output.'
#        sys.exit()
#else:
#        os.mkdir(parameters['output_directory'])

# ----------------------------------------
# 2) LOAD IN NECESSARY FUNCTIONS FROM MODULE FILES
# ----------------------------------------
if parameters['hessian_algorithm'].lower() == 'henm':
        calculate_hessian = importlib.import_module(parameters['hessian_functions_file'].split('.')[0],package=None).perform_henm
#elif parameters['hessian_algorithm'].lower() == 'reach':
#        calculate_hessian = importlib.import_module(parameters['hessian_functions_file'].split('.')[0],package=None).perform_reach

# ----------------------------------------
# MAIN
# ----------------------------------------
if __name__ == '__main__':
	main()

