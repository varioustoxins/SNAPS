import sys
import os
import re
import shutil

from flask import Flask, render_template, jsonify, request
from os import environ
from validation import Validate
from args import Args

mainNAPSfilePath = os.path.dirname(os.path.realpath(__file__)) + '/../python'
sys.path.append(mainNAPSfilePath)
os.chdir(mainNAPSfilePath)

from NAPS import runNAPS
app = Flask(__name__)

@app.route('/run', methods = ['POST'])
def run():
    args = getArgs(request, Args(app.instance_path))
    result = Validate(args)
    if (result[0]):
        try:
            runNAPS(args.argsToList())
            response = createJSONForTable(args.output_file)
        except:
            #log errors
            response = jsonify(status='application_failed')
    else:
        response = result[1]

    if os.path.exists(args.directory) and os.path.isdir(args.directory):
        shutil.rmtree(args.directory)

    return response

def getArgs(request, args):
    #For now, default files are used if files are not provided
    if 'observedShiftsFile' in request.form:
        args.shift_file = '../data/P3a_L273R/naps_shifts.txt'
    else:
        request.files['observedShiftsFile'].save(args.shift_file)

    if 'predictedShiftsFile' in request.form:
        args.pred_file = '../data/P3a_L273R/shiftx2.cs'
    else:
        request.files['predictedShiftsFile'].save(args.pred_file)

    args.shift_type = request.form['shift_type'].strip().lower()
    args.pred_type = request.form['pred_type'].strip().lower()
    return args

def createJSONForTable(output_path):
    with open(output_path) as output_file:
        result = []
        line = output_file.readline()
        headers = re.split(r'\t', line.rstrip('\n'))
        line = output_file.readline()
        while line:
            row = {}
            values = re.split(r'\t', line.rstrip('\n'))
            for i, header in enumerate(headers):
                row[header] = values[i]
            result.append(row)
            line = output_file.readline()
        return jsonify(status='ok', headers=headers, result=result)

@app.route('/')
def index():
    return render_template('index.html')

if __name__ == '__main__':
    HOST = environ.get('SERVER_HOST', 'localhost')
    try:
        PORT = int(environ.get('SERVER_PORT', '5555'))
    except ValueError:
        PORT = 5555
    app.run(HOST, PORT)
