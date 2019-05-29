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
# Configuration file for NAPS
prob_method     pdf     # Method for calculating probability (options are cdf or pdf)
pred_correction False   # Apply a linear correction to the predicted shifts
pred_correction_file      ../config/lin_model_shiftx2.csv    # File containing parameters for linear correction to predicted shift
delta_correlation_mean_file     ../config/d_mean.csv       # File containing mean prediction errors
delta_correlation_cov_file      ../config/d_cov.csv        # File containing covariances between the prediction errors
delta_correlation_mean_corrected_file     ../config/dd_mean.csv       # File containing mean prediction errors, assuming the predictions have been corrected
delta_correlation_cov_corrected_file      ../config/dd_cov.csv        # File containing covariances between the prediction errors, assuming the predictions have been corrected
alt_assignments 0       # Number of alternative assignments to generate
#atom_sd "H:0.1711, N:1.1169, HA:0.1231, C:0.5330, CA:0.4412, CB:0.5163, C_m1:0.5530, CA_m1:0.4412, CB_m1:0.5163"    # Atom standard deviations. Comma separated.
atom_sd "H:0.454, N:2.429, HA:0.227, C:1.030, CA:0.932, CB:1.025, C_m1:1.030, CA_m1:0.932, CB_m1:1.025"    # Atom standard deviations. Comma separated.
plot_strips     True\n"""

    # Write the user-submitted parameters
    s += "pred_offset %s\n" % str(request.form.get("predResOffset"))
    if request.form.get("deltaCorrelation")=="on":
        s += "delta_correlation       True\n"
    else:
        s += "delta_correlation       False\n"
    s += "use_ss_class_info       False\n"
    s += """atom_set      "%s"\n""" % ",".join(request.form.getlist("atomType"))
    #print(s)
    f.write(s)

    f.close()
    