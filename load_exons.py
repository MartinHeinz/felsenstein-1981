import pickle

with open('nonexon_windows_w100_monkeys.pkl', 'rb') as windows:
    alpha_freq_table = pickle.load(windows)

import matplotlib.pyplot as plt

plt.rcdefaults()
import numpy as np

objects = list(map(str, alpha_freq_table.keys()))
y_pos = np.arange(len(objects))
performance = list(alpha_freq_table.values())

fig, ax = plt.subplots()
plt.bar(y_pos, performance, align='center', alpha=0.5, color='b')
plt.xticks(y_pos, objects)
plt.ylabel('Count')
plt.xlabel('Alpha')
plt.title('Non-exon windows only for primates.')
# plt.locator_params(axis='x', tight=True, nbins=4)


plt.show()
