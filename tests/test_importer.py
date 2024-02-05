import pandas as pd
import pytest
from numpy import NaN
from pandas.testing import assert_frame_equal

from SNAPS_assigner import SNAPS_assigner
from SNAPS_importer import SNAPS_importer, SnapsImportException


def test_import_aa_type_info():
    importer = SNAPS_importer()
    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')
    result = importer.import_aa_type_info('test_data/test_aa_info.txt')
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
    result_out = importer.import_aa_type_info('test_data/test_aa_info_out.txt')
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

    expected_out = pd.DataFrame(data_out, index=data_out['SS_name'], columns=['SS_name', 'C', 'CA', 'CB', 'H', 'HA',
                                                                              'N', 'SS_class'])
    expected_out.columns.names = ['Atom_type']
    assert_frame_equal(expected_out, result_out, check_exact=False, rtol=0.01)


def test_headings_aa_info():
    importer = SNAPS_importer()
    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')

    with pytest.raises(SnapsImportException) as e:
        importer.import_aa_type_info('test_data/test_aa_info_bad_header.txt')

    assert 'Unexpected column name(s) [AAA]' in str(e.value)
    assert 'expected column names are:' in str(e.value)
    for name in 'Type SS_name AA'.split():
        assert name in str(e.value)


def test_headings_aa_info_type_column_bad():
    importer = SNAPS_importer()
    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')

    with pytest.raises(SnapsImportException) as e:
        importer.import_aa_type_info('test_data/test_aa_info_bad_type_column.txt')
    assert "Type column row error: 'Type' column rows can only contain 'in' or 'ex'" in str(e.value)


def test_lowercase_headings_bad():
    importer = SNAPS_importer()
    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')

    with pytest.raises(SnapsImportException) as e:
        importer.import_aa_type_info('test_data/test_aa_info_lowercase_column_names.txt')


def test_aa_letters_bad():
    importer = SNAPS_importer()
    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')

    with pytest.raises(SnapsImportException) as e:
        importer.import_aa_type_info('test_data/test_aa_info_letters_bad.txt')
    assert "AA type letters incorrect, Amino Acid letters can only be: 'ACDEFGHIKLMNPQRSTVWY'" in str(e.value)


def test_obs_data_bad():
    importer = SNAPS_importer()
    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')

    with pytest.raises(SnapsImportException) as e:
        importer.import_aa_type_info('test_data/test_aa_info_bad_obs_data.txt')
    assert "Incorrect data given, unexpected spin system in aa types. the input spin systems should be in the chemical"\
           " shift list" in str(e.value)


def test_offset_data():
    importer = SNAPS_importer()
    importer.import_obs_shifts('test_data/nef_resonances_4_test.nef', 'nef')

    with pytest.raises(SnapsImportException) as e:
        importer.import_aa_type_info('test_data/test_aa_info_m1.txt')
    assert "offset can only be 0 or -1" in str(e.value)


# def test_import_P3a_L273R_shifts():
#     importer = SNAPS_importer()
#     importer.import_obs_shifts('test_data/P3a_L273R_10_test.nef', 'nef')
#     result = importer.import_aa_type_info('test_data/sparky_10_test.txt.txt', 'i')

