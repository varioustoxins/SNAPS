from flask import jsonify

class Validate:
    """validate NAPS args"""

    shift_type_error = 'Invalid observed shift type.'
    pred_type_error = 'Invalid predicted shift type.'

    validShiftTypes = [
        'naps',
        'ccpn',
        'sparky',
        'xeasy',
        'nmrpipe',
        'test'
        ]

    validPredTypes = [
        'shiftx2',
        'sparta+'
        ]

    def __new__(self, args):
        errors = []
        if args.shift_type not in self.validShiftTypes:
            errors.append(shift_type_error)
        if args.pred_type not in self.validPredTypes:
            errors.append(pred_type_error)

        isValid = not any(errors)
        return isValid, jsonify(status='validation_failed', errors=errors)