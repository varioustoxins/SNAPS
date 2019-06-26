import uuid
import os
from pathlib import Path
from flask import send_file

class Args:
    """args for NAPS"""

    def __init__(self, instance_path, form):
        self.directory = os.path.join(instance_path, 'tmp_' + str(uuid.uuid4()))
        self.shift_file = os.path.join(self.directory, 'shift.txt')
        self.pred_file = os.path.join(self.directory, 'pred.txt')
        self.config_file = os.path.join(self.directory, 'config.txt')
        self.output_file = os.path.join(self.directory, 'results.txt')
        self.log_file = os.path.join(self.directory, 'log.txt')
        self.plot = {}
        self.plot_file = os.path.join(self.directory, 'plot.html')
        self.shift_type = form['shift_type'].strip().lower()
        self.pred_type = form['pred_type'].strip().lower()
        self.output_shiftlist = os.path.join(self.directory, 'output_shiftlist.txt')
        self.shift_output_type = form['outShiftType'].strip().lower()

    def argsToList(self):
        arg_list = [
            self.shift_file,
            self.pred_file,
            self.output_file,
            '--shift_type', self.shift_type,
            '--pred_type', self.pred_type,
            '--plot_file', self.plot_file,
            '--shift_output_type', self.shift_output_type,
            '--shift_output_file', self.output_shiftlist,
            '-c', self.config_file,
            '-l', self.log_file
        ]
        return arg_list

    def getResults(self):
        with open(self.output_file,mode='r') as f:
            return f.read()
        
    def getShiftlist(self):
        with open(self.output_shiftlist,mode='r') as f:
            return f.read()
        
    def getLog(self):
        with open(self.log_file,mode='r') as f:
            return f.read()

    def getPlot(self):
        hmm = self.plot if self.plot else ''
        return hmm

    def getFiles(self):
        return {
            'results': self.getResults(),
            'shiftlist':self.getShiftlist(),
            'plot': self.getPlot(),
            #'log': self.getLog()
            
        }