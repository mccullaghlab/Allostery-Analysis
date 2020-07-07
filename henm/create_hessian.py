
# ----------------------------------------
# USAGE AND REFERENCING
# ----------------------------------------
#   python3 create_hessian.py create_hessian_config_file_name IO_functions_file

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

# ----------------------------------------
# VARIABLE DECLARATION: 
# ----------------------------------------

config_file = sys.argv[1]
IO_functions_file = sys.argv[2]

config_parser = importlib.import_module(IO_functions_file.split('.py')[0],package=None).create_hessian_config_parser
summary = importlib.import_module(IO_functions_file.split('.py')[0],package=None).create_hessian_summary

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
        if parameters['hessian_algorithm'].lower() == 'henm':   # other hessian model algorithms may be added in the future.
                henm_hessian = calculate_hessian(covariance, average_structure, parameters['output_directory'], guess = guess_hessian, max_iterations = int(parameters['henm_max_iterations']), alpha = float(parameters['henm_alpha']), kBT = 1.9872E-3*float(parameters['temperature']), distance_cutoff = float(parameters['distance_cutoff']), threshold = float(parameters['henm_threshold']))
                np.savetxt(parameters['output_directory'] + parameters['output_hessian_file_name'],henm_hessian)

        if parameters['summary_boolean']:
                summary(parameters['output_directory'] + 'create_hessian.summary',sys.argv,parameters)

# ----------------------------------------
# 1) LOAD IN USER DEFINED PARAMETERS
# ----------------------------------------
parameters = {}
config_parser(config_file,parameters)

# ----------------------------------------
# SETTING UP THE OUTPUT DIRECTORY
# ----------------------------------------
if parameters['output_directory'][-1] != os.sep:
        parameters['output_directory'] += os.sep

# ----------------------------------------
# 2) LOAD IN NECESSARY FUNCTIONS FROM MODULE FILES
# ----------------------------------------
if parameters['hessian_algorithm'].lower() == 'henm':
        calculate_hessian = importlib.import_module(parameters['hessian_functions_file'].split('.')[0],package=None).perform_henm

# ----------------------------------------
# MAIN
# ----------------------------------------
if __name__ == '__main__':
	main()

