# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 12:04:35 2018

@author: Alex
"""

from flask import (
    Flask, Blueprint, flash, g, redirect, render_template, request, session, url_for
)

app = Flask(__name__)
bp = Blueprint('main', __name__, url_prefix='/')
 
@app.route("/")
def index():
    return render_template(
        'index.html',**locals())
    
    
@bp.route('/run', methods=('GET, POST'))
def run():
    print("running")
    file=open("test.txt", "w")
    file.close()
    return render_template(
        'index.html',**locals())
 
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=80)