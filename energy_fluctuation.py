import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

trk = open('data/counter', 'r')
counter = int(trk.read()) - 1

filename = 'data/energy_fluctuation' + str(counter)
data = pd.read_table(filename, sep=r'\s+')

sns.set(context='notebook', style='darkgrid')
sns.set_palette('Paired', color_codes=True)

fig, ax = plt.subplots()
fig.suptitle(r'Fluctuation of energy with respect to dt')

ax.scatter(data['dt'], data['fluctuation'], color='r')
ax.set_xlabel(r"dt $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$")
ax.set_ylabel(r"$\sigma_{energy}$ $\left(\varepsilon\right)$")
# ax.set_xlim(left=10**(-4))
# ax.set_xlim(right=10**(-2))
# ax.set_ylim(top=1)
# ax.set_ylim(bottom=0.1)
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.set_ylim(top=100)
# ax.set_ylim(bottom=0)

fig.set_size_inches(9, 7)
fig.savefig("media/energy_fluctuation" + str(counter) + ".png")
fig.savefig("media/energy_fluctuation" + str(counter) + ".svg")

plt.show()
