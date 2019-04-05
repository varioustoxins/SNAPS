import sys
import os
import re
import base64
import json

from flask import Flask, render_template, jsonify, request, session
from flask_mail import Mail
from os import environ
from validation import Validate
from args import Args
from fileSender import emailFiles
from fileHandler import saveFiles, deleteFiles

mainNAPSfilePath = os.path.dirname(os.path.realpath(__file__)) + '/../python'
sys.path.append(mainNAPSfilePath)
os.chdir(mainNAPSfilePath)

from NAPS import runNAPS
app = Flask(__name__)
app.secret_key = 'napsnapsnapsnaps' #should be changed to an external config value in production

mail_settings = {
    "MAIL_SERVER": '',
    "MAIL_PORT": 0,
    "MAIL_USE_TLS": False,
    "MAIL_USE_SSL": True,
    "MAIL_USERNAME": '',
    "MAIL_PASSWORD": ''
}

app.config.update(mail_settings)
mail = Mail(app)

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
        args.plot = runNAPS(args.argsToList())
        return createJSONForTable(args)
    except Exception as e:
        #log errors
        print("Unexpected error:" + str(e))
        return jsonify(status='application_failed')

def createJSONForTable(args):
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

    return jsonify(status='ok', headers=headers, result=result, files=args.getFiles())

@app.route('/')
def index():
    return render_template('index.html')
    
@app.route('/info')
def info():
    return render_template('info.html')

#@app.route('/email', methods = ['POST'])
#def email():
#    if not app.config.get("MAIL_SERVER"):
#        return jsonify(status='no_email_server', message='The email server is not configured on the views.py page.')
#    return emailFiles(request, mail, app)

if __name__ == '__main__':
    HOST = environ.get('SERVER_HOST', 'localhost')
    try:
        PORT = int(environ.get('SERVER_PORT', '5555'))
    except ValueError:
        PORT = 5555
    app.run(HOST, PORT)
