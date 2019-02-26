import unittest
from validation import Validate
from args import Args

class Validation_test(unittest.TestCase):
    def validate_allInvalid_returnErrors(self):
        args = Args('')
        response = Validate(args)

        self.assertFalse(response[0])
        self.assertIn('errors', response)
        self.assertIn(Validate.shift_type_error, response[1].errors)
        self.assertIn(Validate.pred_type_error, response[1].errors)

    def validate_valid_noErrors(self):
        args = Args('')
        args.shift_type = Validate.validShiftTypes[0]
        args.pred_type = Validate.validPredTypes[0]

        response = Validate(args)

        self.assertTrue(response[0])
        self.assertIn('errors', response)
        self.assertNotIn(Validate.shift_type_error, response[1].errors)
        self.assertNotIn(Validate.pred_type_error, response[1].errors)

if __name__ == '__main__':
    unittest.main()
