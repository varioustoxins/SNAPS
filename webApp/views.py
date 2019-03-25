import sys
import os
import re
import base64
import json

from flask import Flask, render_template, jsonify, request, session, send_file
from os import environ
from validation import Validate
from args import Args
from fileHandler import saveFiles, deleteFiles

mainNAPSfilePath = os.path.dirname(os.path.realpath(__file__)) + '/../python'
sys.path.append(mainNAPSfilePath)
os.chdir(mainNAPSfilePath)

from NAPS import runNAPS
app = Flask(__name__)
app.secret_key = 'napsnapsnapsnaps'

@app.route('/run', methods = ['POST'])
def run():
    args = Args(app.instance_path, request.form)
    saveFiles(request, args)
    validationResult = Validate(args)
    response = run(args) if validationResult.isValid else validationResult.response
    #deleteFiles(args)
    return response

def run(args):
    try:
        args.plot = runNAPS(args.argsToList())
        files = {'results': args.getResultsLocation(), 'plot': args.getPlotLocation()}
        saveSession(files)
        return createJSONForTable(args, files)
    except Exception as e:
        #log errors
        print("Unexpected error:" + str(e))
        return jsonify(status='application_failed')

def createJSONForTable(args, files):
    with open(args.output_file) as output_file:
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

    return jsonify(status='ok', headers=headers, result=result, plot=args.plot, files=files)

@app.route('/')
def index():
    files = getSessionInfo()
    return render_template('index.html', files=files)
    
@app.route('/info')
def info():
    files = getSessionInfo()
    return render_template('info.html', files=files)

@app.route('/download', methods = ['POST'])
def download():
    try:
        return send_file(os.path.join(app.instance_path, request.values['filePath']), as_attachment=True)
    except Exception as e:
        return str(e)

def saveSession(sessionInfo):
    session['userSession'] = sessionInfo

def getSessionInfo():
    if 'userSession' in session:
        return session['userSession']
    else:
        return {}


if __name__ == '__main__':
    HOST = environ.get('SERVER_HOST', 'localhost')
    try:
        PORT = int(environ.get('SERVER_PORT', '5555'))
    except ValueError:
        PORT = 5555
    app.run(HOST, PORT)
