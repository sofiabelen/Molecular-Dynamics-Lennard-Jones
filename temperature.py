import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

trk = open('Data2/counter', 'r')
counter = int(trk.read()) - 1

filename = 'Data2/temperature' + str(counter)
data = pd.read_table(filename, sep=r'\s+')

sns.set(context='notebook', style='darkgrid')
sns.set_palette('Reds', color_codes=True)

fig, ax = plt.subplots()

ax.scatter(data['time'], data['temp'], color='g')

ax.set_xlabel(r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$")
ax.set_ylabel(r"Темература $\left(\varepsilon\right)$")
ax.set_title(r"Зависимость температуры от времени")

fig.set_size_inches(8, 6)
plt.savefig("Data2/temperature" + str(counter)+".png")
