import sys
import os
import re

from flask import Flask, render_template, jsonify, request
from os import environ
from validation import Validate
from args import Args
from fileHandler import saveFiles, deleteFiles

mainNAPSfilePath = os.path.dirname(os.path.realpath(__file__)) + '/../python'
sys.path.append(mainNAPSfilePath)
os.chdir(mainNAPSfilePath)

from NAPS import runNAPS
app = Flask(__name__)

@app.route('/run', methods = ['POST'])
def run():
    args = Args(app.instance_path, request.form)
    saveFiles(request, args)
    validationResult = Validate(args)
    response = run(args) if validationResult.isValid else validationResult.response
    deleteFiles(args)
    return response

def run(args):
    try:
        runNAPS(args.argsToList())
        return createJSONForTable(args.output_file)
    except:
        #log errors
        return jsonify(status='application_failed')

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
