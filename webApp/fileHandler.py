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

def deleteFiles(args):
    if os.path.exists(args.directory) and os.path.isdir(args.directory):
        shutil.rmtree(args.directory)
