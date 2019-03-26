import uuid
import os

class Args:
    """args for NAPS"""

    def __init__(self, instance_path, form):
        self.directory = os.path.join(instance_path, 'tmp_' + str(uuid.uuid4()))
        self.shift_file = os.path.join(self.directory, 'shift.txt')
        self.pred_file = os.path.join(self.directory, 'pred.txt')
        self.output_file = os.path.join(self.directory, 'results.txt')
        self.plot = {}
        self.plot_file = os.path.join(self.directory, 'plot.html')
        self.shift_type = form['shift_type'].strip().lower()
        self.pred_type = form['pred_type'].strip().lower()

    def argsToList(self):
        return [
            self.shift_file,
            self.pred_file,
            self.output_file,
            '--shift_type', self.shift_type,
            '--pred_type', self.pred_type,
            '--plot_file', self.plot_file,
            '-c', '../config/config.txt',
            #'-l', '../output/test.log'
        ]

    def getResultsLocation(self):
        return self.output_file

    def getPlotLocation(self):
        return self.plot_file if self.plot else ''