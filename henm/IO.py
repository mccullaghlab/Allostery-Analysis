
# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------
def create_hessian_config_parser(config_file,parameters):	# Function to take config file and create/fill the parameter dictionary 
        """ Function to take config file and create/fill the parameter dictionary (created before function call). 
        
        Usage: 
            parameters = {}     # initialize the dictionary to be filled with keys and values
            trajectory_analysis_config_parser(config_file,parameters)

        Arguments:
            config_file: string object that corresponds to the local or global position of the config file to be used for this analysis.

        """
        necessary_parameters = ['output_directory','hessian_algorithm','system_covariance_file','system_average_structure_file','hessian_functions_file']

        all_parameters = ['output_directory','hessian_algorithm','system_covariance_file','system_average_structure_file','hessian_functions_file','summary_boolean','plotting_boolean','initial_guess_hessian','henm_alpha','henm_max_iterations','henm_threshold','temperature']
        
        for i in range(len(necessary_parameters)):
            parameters[necessary_parameters[i]] = ''
        
        # SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['summary_boolean'] = False 
        parameters['plotting_boolean'] = False 
        parameters['initial_guess_hessian'] = None
        parameters['henm_alpha'] = 1E-2
        parameters['henm_max_iterations'] = 100
        parameters['henm_threshold'] = 1E-4
        parameters['temperature'] = 303.
        #parameters[''] = 
        
        # GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
        with open(config_file) as f:
            exec(compile(f.read(),config_file,'exec'),parameters)
        
        for key, value in list(parameters.items()):
            if value == '':
                print('%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key))
                sys.exit()

def create_hessian_summary(summary_file_name,arguments,parameters):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            create_hessian_analysis_summary(summary_file_name,arguments,parameters)

        Arguments:
            summary_file_name: string object of the file name to be written that holds the summary information.
            parameters: dictionary object filled with the parameters used in the analysis.

        """
        with open(summary_file_name,'w') as f:
            f.write('To recreate this analysis, run this line:\n')
            for i in range(len(arguments)):
                f.write('%s ' %(arguments[i]))
            f.write('\n\n')
            f.write('Parameters used:\n')
            for key, value in list(parameters.items()):
                if key == '__builtins__':
                    continue
                if type(value) == int or type(value) == float:
                    f.write("%s = %s\n" %(key,value))
                else:
                    f.write("%s = '%s'\n" %(key,value))

