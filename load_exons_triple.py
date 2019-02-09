import pickle

with open('exon_windows_w25.pkl', 'rb') as windows:
    alpha_freq_table25 = pickle.load(windows)

with open('exon_windows_w200.pkl', 'rb') as windows:
    alpha_freq_table200 = pickle.load(windows)

with open('exon_windows_w500.pkl', 'rb') as windows:
    alpha_freq_table500 = pickle.load(windows)

import matplotlib.pyplot as plt

plt.rcdefaults()
import numpy as np
from matplotlib.ticker import FuncFormatter

formatter = FuncFormatter(lambda y, pos: "%d%%" % (y))

y = list(map(str, alpha_freq_table25.keys()))
y_pos = np.arange(len(y))
x1 = list(alpha_freq_table25.values())
x2 = list(alpha_freq_table200.values())
x3 = list(alpha_freq_table500.values())

x1_count = sum(list(alpha_freq_table25.values()))
x2_count = sum(list(alpha_freq_table200.values()))
x3_count = sum(list(alpha_freq_table500.values()))

x1a = list(map(lambda x: float(x)/x1_count*100, x1))
x2a = list(map(lambda x: float(x)/x2_count*100, x2))
x3a = list(map(lambda x: float(x)/x3_count*100, x3))

# ax = plt.subplot(111)
fig, ax = plt.subplots()
plt.bar(y_pos + 0.00, x1a, color='b', width=0.25, label="w = 25")
plt.bar(y_pos + 0.25, x2a, color='g', width=0.25, label="w = 200")
plt.bar(y_pos + 0.50, x3a, color='r', width=0.25, label="w = 500")
plt.xticks(y_pos, y)
ax.yaxis.set_major_formatter(formatter)
plt.ylabel('Percentage')
plt.xlabel('Alpha')
plt.title('Percentage of exon windows per Alpha')
# plt.locator_params(axis='x', tight=True, nbins=4)
plt.legend()

plt.show()
