import unittest
import json
from .validation import Validate
from .args import Args
from flask import Flask
app = Flask(__name__)

class Tests_Validation(unittest.TestCase):
    def setUp(self):
        self.app_context = app.app_context()
        self.app_context.push()

    def test_validate_allInvalid_returnAllErrors(self):
        dummyForm = {
            'shift_type':'invalid',
            'pred_type':'invalid',
            }
        args = Args('', dummyForm)
        validationResult = Validate(args)

        response = json.loads((validationResult.response.data).decode('utf8'))

        self.assertFalse(validationResult.isValid)
        self.assertIn('errors', response)
        self.assertIn(Validate.shift_type_error, response['errors'])
        self.assertIn(Validate.pred_type_error, response['errors'])

    def test_validate_valid_noErrors(self):
        dummyForm = {
            'shift_type':Validate.validShiftTypes[0],
            'pred_type':Validate.validPredTypes[0],
            }
        args = Args('', dummyForm)
        validationResult = Validate(args)

        response = json.loads((validationResult.response.data).decode('utf8'))

        self.assertTrue(validationResult.isValid)
        self.assertIn('errors', response)
        self.assertNotIn(Validate.shift_type_error, response['errors'])
        self.assertNotIn(Validate.pred_type_error, response['errors'])

if __name__ == '__main__':
    unittest.main()
