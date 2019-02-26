import time
import os

class Args:
    """args for NAPS"""

    def __init__(self, instance_path):
        """ Create a new point at the origin """
        self.directory = os.path.join(instance_path, 'tmp_' + str(time.time()))
        if instance_path:
            os.makedirs(self.directory, exist_ok=True)
        self.shift_file = os.path.join(self.directory, 'shift.txt')
        self.pred_file = os.path.join(self.directory, 'pred.txt')
        self.output_file = os.path.join(self.directory, 'output.txt')
        self.shift_type = ''
        self.pred_type = ''

    def argsToList(self):
        return [
            self.shift_file,
            self.pred_file,
            self.output_file,
            '--shift_type', self.shift_type,
            '--pred_type', self.pred_type,
            '-c', '../config/config.txt',
            #'-l', '../output/test.log'
        ]