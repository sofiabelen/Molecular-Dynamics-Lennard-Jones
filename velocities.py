import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

trk = open('Data2/counter', 'r')
counter = int(trk.read()) - 1

filename = 'Data2/velocities' + str(counter)

data = pd.read_table(filename, sep=r'\s+')

data = data * data

sns.set(context='notebook', style='darkgrid')
sns.set_palette(sns.color_palette("Blues", 3))
# sns.set_palette("Paired")

fig, ax = plt.subplots()

ax = sns.distplot(data['0'], label='x-компонент')
ax = sns.distplot(data['1'], label='y-компонент')
ax = sns.distplot(data['2'], label='z-компонент')

ax.set_xlim(left=0)
ax.set_yscale('log')

ax.set_title(r'Распределение скоростей по проекциям')
ax.set_xlabel(r"$v_i^2$")
ax.set_ylabel(r"$\ln(n_v)$")
ax.legend()

fig.set_size_inches(8, 6)
fig.savefig("Images/velocities" + str(counter) + ".png")
fig.savefig("Images/velocities" + str(counter) + ".svg")
