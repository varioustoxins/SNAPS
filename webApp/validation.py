from flask import jsonify

class Validate:
    """validate NAPS args"""

    shift_type_error = 'Invalid observed shift type.'
    pred_type_error = 'Invalid predicted shift type.'

    validShiftTypes = [
        'snaps',
        'ccpn',
        'sparky',
        'xeasy',
        'nmrpipe',
        'mars',
        'test'
        ]

    validPredTypes = [
        'shiftx2',
        'sparta+'
        ]

    def __new__(self, args):
        errors = []
        if args.shift_type not in self.validShiftTypes:
            errors.append(self.shift_type_error)
        if args.pred_type not in self.validPredTypes:
            errors.append(self.pred_type_error)

        isValid = not any(errors)
        return ValidationResult(isValid, jsonify(status='validation_failed', errors=errors))

class ValidationResult:
    def __init__(self, isValid, response):
        self.isValid = isValid
        self.response = response