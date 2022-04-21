import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

trk = open('data/counter', 'r')
counter = int(trk.read()) - 1

filename = 'data/temperature' + str(counter)
data = pd.read_table(filename, sep=r'\s+')

sns.set(context='notebook', style='darkgrid')
sns.set_palette('Paired', color_codes=True)

fig, ax = plt.subplots()

ax.scatter(data['time'], data['temp'], color='g')

ax.set_xlabel(r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$")
ax.set_ylabel(r"Темература $\left(\varepsilon\right)$")
ax.set_title(r"Зависимость температуры от времени")
ax.set_xlim(left=0)

fig.set_size_inches(8, 6)
plt.savefig("media/temperature" + str(counter)+".png")
plt.savefig("media/temperature" + str(counter)+".svg")
