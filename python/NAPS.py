#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Main NAPS script for assigning an observed shift list based on predicted shifts

@author: aph516
"""
def runNAPS(args):
    #import pandas as pd
    from NAPS_importer import NAPS_importer
    from NAPS_assigner import NAPS_assigner
    import argparse
    #from pathlib import Path
    import logging

    #### Command line arguments
    parser = argparse.ArgumentParser(description="NAPS (NMR Assignments from Predicted Shifts)")
    parser.add_argument("shift_file", 
                        help="A table of observed chemical shifts.")
    parser.add_argument("pred_file",
                        help="A table of predicted chemical shifts.")
    parser.add_argument("output_file",
                        help="The file results will be written to.")

    parser.add_argument("--shift_type", 
                        choices=["naps", "ccpn", "sparky", 
                                 "xeasy", "nmrpipe", "test"], 
                        default="naps", 
                        help="The format of the observed shift file.")
    parser.add_argument("--pred_type", 
                        choices=["shiftx2", "sparta+"],
                        default="shiftx2", 
                        help="The format of the predicted shifts")
    parser.add_argument("-c", "--config_file", 
                        default="/Users/aph516/GitHub/NAPS/python/config.txt",
                        help="A file containing parameters for the analysis.")
    parser.add_argument("-l", "--log_file", default=None,
                        help="A file logging information will be written to.")
    
    parser.add_argument("--shift_output_file", default=None,
                        help="""The file the assigned shiftlist will be written to.""")
    parser.add_argument("--shift_output_type", default="sparky",
                        choices=["sparky","xeasy","nmrpipe"],
                        help="One or more output formats for chemical shift export")
    parser.add_argument("--shift_output_confidence", nargs="*",
                        choices=["High","Medium","Low","Unreliable","Undefined"],
                        default=["High","Medium","Low","Unreliable","Undefined"],
                        help="""Limits the shiftlist output to assignments with 
                        particular confidence levels. More than one level is allowed""")
        
    parser.add_argument("--strip_plot_file", 
                        default=None,
                        help="A filename for an output strip plot.")
    parser.add_argument("--hsqc_plot_file", 
                        default=None,
                        help="A filename for an output HSQC plot.")
    
    #parser.add_argument("--delta_correlation", action="store_true", 
    #                    help="If set, account for correlations between prediction errors of different atom types")
    parser.add_argument("-a", "--alt_assignments", default=-1, type=int,
                        help="The number of alternative assignments to generate, "+
                        "in addition to the highest ranked.")
    parser.add_argument("--test_aa_classes", default=None, 
                        help="""For test data only. 
                        A string containing a comma-separated list of the amino acid 
                        classes for the i residue, a semicolon, then a list of AA 
                        classes for the i-1 residue. No spaces. 
                        eg. "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML" for 
                        a sequential HADAMAC """)
    parser.add_argument("--iterated", action="store_true", 
                        help="If set, use iterated procedure to resolve bad sequential links")
    parser.add_argument("--consistent", action="store_true", 
                        help="If set, enforce consistent links for high confidence assignments and iterate.")
    
    
    if True:
        args = parser.parse_args(args)
    else:   # For testing
        args = parser.parse_args(("../data/P3a_L273R/naps_shifts.txt",
                                  "../data/P3a_L273R/shiftx2.cs",
                                  "../output/test.txt",
                                  "--shift_type","naps",
                                  "--pred_type","shiftx2",
                                  "-c","../config/config.txt",
                                  "-l","../output/test.log"))    

    # Set up logging
    if args.log_file is not None:
        # Create a logger
        logger = logging.getLogger("SNAPS")
        logger.setLevel(logging.DEBUG)
        # Create a log handler that writes to a specific file.
        # In principle you could have multiple handlers, but here I just have one.
        # Need to explicitly define a handler so it can be explicitly closed 
        # once the analysis is complete.
        log_handler = logging.FileHandler(args.log_file, mode='w')
        log_handler.setLevel(logging.DEBUG)
        log_handler.setFormatter(logging.Formatter("%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s"))
        logger.addHandler(log_handler)

    else:
        logging.basicConfig(level=logging.ERROR)

    #%%
    #### Set up the NAPS_assigner object
    a = NAPS_assigner()

    # Import config file
    a.read_config_file(args.config_file)
    logger.info("Read in configuration from %s.", args.config_file)

    # Account for any command line arguments that overide config file
    if args.alt_assignments>=0:
        a.pars["alt_assignments"] = args.alt_assignments

    # Import observed and predicted shifts
    importer = NAPS_importer()

    if args.shift_type=="test":
        if args.test_aa_classes is None:
            importer.import_testset_shifts(args.shift_file)
        else:
            AA_class, AA_class_m1 = args.test_aa_classes.split(";")
            importer.import_testset_shifts(args.shift_file,
                                           SS_class=AA_class.split(","),
                                           SS_class_m1=AA_class_m1.split(","))
    else:
        importer.import_obs_shifts(args.shift_file, args.shift_type, SS_num=False)
        
    a.obs = importer.obs
    logger.info("Read in %d spin systems from %s.", 
                 len(a.obs["SS_name"]), args.shift_file)


    a.import_pred_shifts(args.pred_file, args.pred_type)
    logger.info("Read in %d predicted residues from %s.", 
                 len(a.preds["Res_name"]), args.pred_file)

    #### Do the analysis
#    if args.shift_type=="test":
#        # Limit observations and predictions to their common range
##        a.obs = a.obs[a.obs["Res_N"]>=a.preds["Res_N"].min()]
##        a.obs = a.obs[a.obs["Res_N"]<=a.preds["Res_N"].max()]
##        a.preds = a.preds[a.preds["Res_N"]>=a.obs["Res_N"].min()]
##        a.preds = a.preds[a.preds["Res_N"]<=a.obs["Res_N"].max()]
#        # Delete any unmatchable observations or predictions
#        a.obs = a.obs[a.obs["SS_name"].isin(a.preds["Res_name"])]
#        a.preds = a.preds[a.preds["Res_name"].isin(a.obs["SS_name"])]
        
    a.add_dummy_rows()
    a.calc_log_prob_matrix2(sf=1, verbose=False)
    logger.info("Calculated log probability matrix (%dx%d).", 
                 a.log_prob_matrix.shape[0], a.log_prob_matrix.shape[1])
    a.calc_mismatch_matrix()
    matching = a.find_best_assignments()
    a.make_assign_df(matching, set_assign_df=True)
    logger.info("Calculated best assignment.")
    a.add_consistency_info(threshold=a.pars["seq_link_threshold"])
    logger.info("Checked assignment consistency.")
    
    #### Output the results
    if args.iterated:
        matching2 = a.find_consistent_assignments(10)
        a.make_assign_df(matching2, set_assign_df=True)
        a.add_consistency_info(threshold=a.pars["seq_link_threshold"])
        a.assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
                           index=False)
    elif args.consistent:
        a.assign_df = a.find_consistent_assignments2()
        a.assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
                           index=False)
    elif a.pars["alt_assignments"]>0:
        a.find_alt_assignments(N=a.pars["alt_assignments"], verbose=False, 
                                by_ss=True)
        logger.info("Calculated the %d next best assignments for each spin system", 
                     a.pars["alt_assignments"])
        a.alt_assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
                               index=False)
    else:
        a.assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
                           index=False)
    
    logger.info("Wrote results to %s", args.output_file)
    
    #### Write chemical shift lists
    if args.shift_output_file is not None:
        a.output_shiftlist(args.shift_output_file, args.shift_output_type, 
                           confidence_list=args.shift_output_confidence)
    
    #### Make some plots
    plots = []
    if args.hsqc_plot_file is not None:
        hsqc_plot = a.plot_hsqc_bokeh(args.hsqc_plot_file, "html")
        logger.info("Wrote HSQC plot to %s", args.hsqc_plot_file)
        plots += [hsqc_plot]
        
    if args.strip_plot_file is not None:
        strip_plot = a.plot_strips_bokeh(args.strip_plot_file, "html")
        logger.info("Wrote strip plot to %s", args.strip_plot_file)
        plots += [strip_plot]
    
    # Close the log file
    logger.handlers[0].close()
    logger.removeHandler(logger.handlers[0])

    return(plots)


#%% Run the actual script
if __name__ == '__main__':
    import sys

    runNAPS(sys.argv[1:])
