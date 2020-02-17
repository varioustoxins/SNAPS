# Flask web app for SNAPS

from flask import Flask, render_template, jsonify, request
#from SNAPS import SNAPS_compute
import sys, os

#main_SNAPS_file_path = os.path.dirname(os.path.realpath(__file__)) + '/../python'
#sys.path.append(main_SNAPS_file_path)
#os.chdir(main_SNAPS_file_path)

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/info")
def info():
    return render_template("info.html")

@app.route("/howto")
def howto():
    return render_template("howto.html")

if __name__ == '__main__':
    app.run(debug=True)