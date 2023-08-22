import pandas as pd
import pytest
from numpy import NaN
from pandas.testing import assert_frame_equal
from SNAPS_importer import SNAPS_importer, SnapsImportException


def test_import_aa_type_info():
    importer = SNAPS_importer()

    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')
    result = importer.import_aa_type_info('test_data/test_aa_info.txt', 'i')

    data = {
        'Atom_type': ['236Pro', '237Ala', '238Met', '239Thr'],
        'SS_name': ['236Pro', '237Ala', '238Met', '239Thr'],
        'C': [176.44282, 177.76473, 175.89544, 173.0364],
        'CA': [62.82192, 52.38738, 55.48441, 62.02754],
        'CB': [32.26901, 19.14546, 32.90399, 69.45566],
        'H': [NaN, 8.54894, 8.37917, 7.88704],
        'HA': [4.43775, 4.27742, 4.41167, 3.88160],
        'N': [NaN, 124.63832, 119.82062, 115.46265],
        'SS_class': ['ACDEFGHIKLMNPQRSTVWY', 'AVI', 'AVI', 'ACDEFGHIKLMNPQRSTVWY']
    }
    expected = pd.DataFrame(data, index=data['SS_name'], columns=['SS_name', 'C', 'CA', 'CB', 'H', 'HA', 'N',
                                                                  'SS_class'])

    expected.columns.names = ['Atom_type']
    assert_frame_equal(expected, result, check_exact=False, rtol=0.01)


def test_import_aa_type_info_out():
    importer = SNAPS_importer()

    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')
    result_out = importer.import_aa_type_info('test_data/test_aa_info_out.txt', 'i')

    data_out = {
        'Atom_type': ['236Pro', '237Ala', '238Met', '239Thr'],
        'SS_name': ['236Pro', '237Ala', '238Met', '239Thr'],
        'C': [176.44282, 177.76473, 175.89544, 173.0364],
        'CA': [62.82192, 52.38738, 55.48441, 62.02754],
        'CB': [32.26901, 19.14546, 32.90399, 69.45566],
        'H': [NaN, 8.54894, 8.37917, 7.88704],
        'HA': [4.43775, 4.27742, 4.41167, 3.88160],
        'N': [NaN, 124.63832, 119.82062, 115.46265],
        'SS_class': ['ACDEFGHIKLMNPQRSTVWY', 'CDEFGHKLMNPQRSTWY', 'CDEFGHKLMNPQRSTWY', 'ACDEFGHIKLMNPQRSTVWY']
    }
    expected_out = pd.DataFrame(data_out, index=data_out['SS_name'], columns=['SS_name', 'C', 'CA', 'CB', 'H', 'HA'
                                                                                                    , 'N', 'SS_class'])
    expected_out.columns.names = ['Atom_type']

    assert_frame_equal(expected_out, result_out, check_exact=False, rtol=0.01)



