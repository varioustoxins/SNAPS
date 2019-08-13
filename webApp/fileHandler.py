import os
import shutil

def saveFiles(request, args):
    os.makedirs(args.directory, exist_ok=True)

    #For now, default files are used if files are not provided
    if 'observedShiftsFile' in request.form:
        args.shift_file = '../data/testset/simplified_BMRB/4834.txt'
        args.shift_type = 'test'
    else:
        request.files['observedShiftsFile'].save(args.shift_file)

    if 'predictedShiftsFile' in request.form:
        args.pred_file = '../data/testset/shiftx2_results/A003_1LM4B.cs'
        args.pred_type = 'shiftx2'
    else:
        request.files['predictedShiftsFile'].save(args.pred_file)
    
    #args.config_file = '../config/config_plot.txt'
    makeConfigFile(request, args)

def deleteFiles(args):
    if os.path.exists(args.directory) and os.path.isdir(args.directory):
        shutil.rmtree(args.directory)

def makeConfigFile(request, args):
    """Create a configuration file using user-submitted options"""
    
    f = open(args.config_file, "w+")
    
    # Write the constant parameters to the file
    s = """
# SNAPS configuration file
# Follows the YAML format
# Use spaces for indentation, not tabs!     
atom_sd: # Atom standard deviations. Don't change indentation.
    H: 0.454 
    N: 2.429 
    HA: 0.227 
    C: 1.030 
    CA: 0.932 
    CB: 1.025
    C_m1: 1.030 
    CA_m1: 0.932 
    CB_m1: 1.025
iterate_until_consistent:       False   # If True, iteratively enforce consistent links for High and Medium confidence assignments
delta_correlation:       True   # Account for correlations in prediction errors
delta_correlation_mean_file:     ../config/d_mean.csv       # File containing mean prediction errors
delta_correlation_cov_file:      ../config/d_cov.csv        # File containing covariances between the prediction errors
pred_correction:        False   # Apply a linear correction to the predicted shifts
pred_correction_file:      ../config/lin_model_shiftx2.csv    # File containing parameters for linear correction to predicted shift
delta_correlation_mean_corrected_file:     ../config/dd_mean.csv       # File containing mean prediction errors, assuming the predictions have been corrected
delta_correlation_cov_corrected_file:      ../config/dd_cov.csv        # File containing covariances between the prediction errors, assuming the predictions have been corrected
"""

    # Write the user-submitted parameters
    if request.form.get("deltaCorrelation")=="on":
        s += "delta_correlation:       True\n"
    else:
        s += "delta_correlation:       False\n"
    s += "atom_set: \n    - %s\n" % "\n    - ".join(request.form.getlist("atomType"))
    s += "seq_link_threshold:    %s\n" % str(request.form.get("seqLinkThreshold"))
    
    #print(s)
    f.write(s)
    f.close()
    