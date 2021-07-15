from NorthNet.Classes import DataReport
import numpy as np


'''
A test for correct file and data parsing using the DataReport object.

prints out if the test is passed or not
'''

expected_report = DataReport()
expected_report.filename = 'Test_data_report.csv'
expected_report.experiment_code = 'test_data'
expected_report.conditions = {'reactor_volume/ uL': 411.0,
                            'flow_profile_time/ s': np.array([0., 2., 4.]),
                            'Residence time/ s': np.array([120., 120., 120.]),
                            'O=C(CO)CO/ M': 0.4,
                            'C=O/ M': 1.6,
                            '[OH-]/ M': 0.12,
                            'O/ M': 0.0,
                            'O=C(CO)CO_flow_rate/ ul/h': np.array([
                                                                    1541.25,
                                                                    1541.25,
                                                                    1541.25
                                                        ]),
                            'C=O_flow_rate/ ul/h': np.array([
                                                                770.625,
                                                                770.625,
                                                                770.625
                                                    ]),
                            '[OH-]_flow_rate/ ul/h': np.array([
                                                                3082.5,
                                                                3082.5,
                                                                3082.5
                                                    ]),
                            'CaCl2_flow_rate/ ul/h': np.array([
                                                                3082.5,
                                                                3082.5,
                                                                3082.5
                                                    ]),
                            'O_flow_rate/ ul/h': np.array([
                                                            3853.125,
                                                            3853.125,
                                                            3853.125
                                                ])
}

expected_report.analysis_details = {'Calibrations_bounds_date': ['2019_02_18'],
                            'Chromatography_method': ['GCMS', 'Method 5'],
                            'Derivatisation_method': [
                                    '75 uL EtONH2.HCl; 30 mins @ 70 oC; 25 uL N'
                                    ]
}
expected_report.series_values = np.array([
                                            9130.25,  9191.75,  9253.25,
                                            9314.75,  9376.25,  9437.75,
                                            9499.25,  9560.75, 9622.25,
                                            9683.75,  9745.25,  9806.75,
                                            9929.75, 9991.25, 10052.75,
                                            10114.25, 10175.75, 10298.75
                                ])
expected_report.series_unit = 'time/ s'
expected_report.data = {
        'OC[C@@H](O)[C@H](O)CO/ M': np.array([
            0.00083438, 0.00091305, 0.00083907, 0.0007385 , 0.00087215,
            0.00082483, 0.00089383, 0.00082638, 0.00112772, 0.00075661,
            0.00088973, 0.00085594, 0.00076052, 0.00084096, 0.0007136 ,
            0.000885  , 0.00078099, 0.00084392
                                    ]),
        'O=C[C@@H](O)[C@H](O)CO/ M': np.array([
            0.00248029, 0.00278809, 0.00279866, 0.00410499, 0.00354283,
            0.00222424, 0.00209485, 0.00353792, 0.00368276, 0.00245254,
            0.00335681, 0.00242988, 0.00259372, 0.00278114, 0.00322447,
            0.00201809, 0.00164878, 0.00175474
                                    ]),
        'O=CC(O)(CO)CO/ M': np.array([
            0.01779828, 0.01761199, 0.01836101, 0.01718358, 0.01778736,
            0.01657012, 0.01700838, 0.01733549, 0.01723092, 0.01716058,
            0.01684084, 0.01828846, 0.01815032, 0.01771439, 0.01836283,
            0.01716965, 0.01668995, 0.01885896
                                    ])
}

expected_report.errors = {
        'OC[C@@H](O)[C@H](O)CO/ M': np.array([
            10.00091305, 10.00083907, 10.0007385 , 10.00087215, 10.00082483,
            10.00089383, 10.00082638, 10.00112772, 10.00075661, 10.00088973,
            10.00085594, 10.00076052, 10.00084096, 10.0007136 , 10.000885  ,
            10.00078099, 10.00084392, 10.
                                    ]),
        'O=C[C@@H](O)[C@H](O)CO/ M': np.array([
            10.00278809, 10.00279866, 10.00410499, 10.00354283, 10.00222424,
            10.00209485, 10.00353792, 10.00368276, 10.00245254, 10.00335681,
            10.00242988, 10.00259373, 10.00278114, 10.00322447, 10.00201809,
            10.00164878, 10.00175474, 10.
                                    ]),
       'O=CC(O)(CO)CO/ M': np.array([
            10.01761199, 10.01836101, 10.01718358, 10.01778736, 10.01657012,
            10.01700838, 10.01733549, 10.01723092, 10.01716058, 10.01684084,
            10.01828846, 10.01815032, 10.01771439, 10.01836283, 10.01716965,
            10.01668995, 10.01885896, 10.
                            ])
}

# Reading from a file
report = DataReport(file = 'tests/Test_data_report.csv')

checks = {}
checks['filename'] = report.filename == expected_report.filename
checks['experiment code'] = report.experiment_code == expected_report.experiment_code

for c in expected_report.conditions:
    # contents are either float or numpy array
    name = 'conditions->{}'.format(c)
    if c not in report.conditions:
        checks[name] = False
    elif type(report.conditions[c]) == float:
        checks[name] = report.conditions[c] == expected_report.conditions[c]
    else:
        checks[name] = np.allclose(
                        report.conditions[c], expected_report.conditions[c],
                        atol = 1e-8
                        )

for a in expected_report.analysis_details:
    # all values are lists containing strings
    name = 'analysis_details->{}'.format(a)
    if a not in report.analysis_details:
        checks[name] = False
    else:
        checks[name] = report.analysis_details[a] == expected_report.analysis_details[a]


checks['series_values'] = np.allclose(
                        report.series_values, expected_report.series_values,
                        atol = 1e-8
                        )

checks['series_unit'] = report.series_unit == expected_report.series_unit

for d in expected_report.data:
    name = 'data->{}'.format(d)
    if d not in report.data:
        checks[name] = False
    else:
        checks[name] = np.allclose(report.data[d], expected_report.data[d], atol = 1e-8)

for e in expected_report.errors:
    name = 'errors->{}'.format(e)
    if e not in report.errors:
        checks[name] = False
    else:
        checks[name] = np.allclose(
                                report.errors[d], expected_report.errors[d],
                                atol = 1e-8
                        )


print(all(list(checks.values())))
