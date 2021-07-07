from NorthNet.Classes import DataReport
import csv


'''
Prints True if the output file writing is identical
to the input text file.
'''
# the output of DataReport.to_string()
# should be the same as read from the raw file
authentic_lines = []
fname = 'tests/Test_data_report.csv'
with open(fname, 'r') as f:
    for line in f:
        ins = line.strip('\n').split(',')
        authentic_lines.append(','.join([x for x in ins if x != '']))

bench_text = "\n".join(authentic_lines)

# Reading from a file
report = DataReport(file = 'tests/Test_data_report.csv')


test_text = report.to_string()


print(bench_text==test_text)
