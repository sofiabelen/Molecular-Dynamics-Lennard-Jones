import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

trk = open('Data2/counter', 'r')
counter = int(trk.read()) - 1

filename = 'Data2/energy' + str(counter)
data = pd.read_table(filename, sep=r'\s+')

sns.set(context='notebook', style='darkgrid')
sns.set_palette('Paired', color_codes=True)

fig, ax = plt.subplots()
ax.scatter(data['time'], data['kinetic'], color='r',\
        label='Средняя кинетическая энергия')
ax.scatter(data['time'], data['potential'], color='g',\
        label='Средняя потенциальная энергия')
ax.scatter(data['time'], data['kinetic'] + data['potential'],\
        color='b', label='Средняя полная энергия')

ax.set_xlabel(r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$")
ax.set_ylabel(r"Энергия $\left(\varepsilon\right)$")
ax.set_title(r"Зависимость энергии от времени")
ax.legend(loc='best')
ax.set_xlim(left=0)

fig.set_size_inches(8, 6)
fig.savefig("Images/energy" + str(counter) + ".png")
fig.savefig("Images/energy" + str(counter) + ".svg")
