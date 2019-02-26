import unittest, os
from webApp import fileHandler
from .args import Args
from flask import Flask
app = Flask(__name__)

class Tests_Files(unittest.TestCase):
    def setUp(self):
        self.app_context = app.app_context()
        self.app_context.push()

    def test_files_createMultipleDirectoriesAndDeleteOne_oneDirectoryRemains(self):
        dummyForm = {'shift_type':'',
         'pred_type':''}
        args1 = Args(os.getcwd(), dummyForm)
        args2 = Args(os.getcwd(), dummyForm)
        try:
            try:
                fileHandler.saveFiles(None, args1)
            except AttributeError:
                self.assertTrue(os.path.exists(args1.directory))

            try:
                fileHandler.saveFiles(None, args2)
            except AttributeError:
                self.assertTrue(os.path.exists(args2.directory))

            fileHandler.deleteFiles(args1)
            self.assertFalse(os.path.exists(args1.directory))
            self.assertTrue(os.path.exists(args2.directory))
            fileHandler.deleteFiles(args2)
            self.assertFalse(os.path.exists(args2.directory))
        finally:
            fileHandler.deleteFiles(args1)
            fileHandler.deleteFiles(args2)

if __name__ == '__main__':
    unittest.main()