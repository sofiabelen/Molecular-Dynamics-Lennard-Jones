import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit

trk = open('Data2/counter', 'r')
counter = int(trk.read()) - 1
filename = 'Data2/positions' + str(counter)

param_name = 'Data2/parameters'  + str(counter)

param = np.loadtxt(param_name)
data = np.loadtxt(filename)

n_steps = int(param[1] - param[2])
n_part = int(param[3])
dim = int(param[4])
dt = int(param[5])

data = data.reshape((n_steps, n_part, dim))

t_range = int((n_steps-1) / 2)

msd = np.zeros(t_range - 1)

for dt in range(1, t_range):
    for d in range(0, dim):
        for i in range(0, n_part):
            for t0 in range(0, t_range):
                msd[dt - 1] +=\
                        (data[t0 + dt][i][d] - data[t0][i][d])**2
    msd[dt - 1] /= float(n_part * t_range * dim)


def line(x, a, b):
    return a*x + b

def parabola(x, a, b, c):
    return a*x**2 + b*x + c

parabola_end = int(t_range / 10)
line_start = int(t_range * 0.7)

popt_parabola, pcov_parabola = curve_fit(f=parabola,\
        xdata=np.arange(0, parabola_end) * dt,\
        ydata=msd[0:parabola_end])

popt_line, pcov_line = curve_fit(f=line,\
        xdata=np.arange(line_start, t_range-1) * dt,\
        ydata=msd[line_start:t_range-1])


sns.set(context='notebook', style='darkgrid')
sns.set_palette('colorblind', color_codes=True)

fig, ax = plt.subplots()
ax.scatter(np.arange(0, t_range - 1) * dt, msd, color='c')

ax.plot(np.arange(0, t_range - 1) * dt,\
        parabola(np.arange(0, t_range - 1) * dt,\
        *popt_parabola), label='Параболическая асимтота',\
        color='b')

ax.plot(np.arange(0, t_range - 1) * dt,\
        line(np.arange(0, t_range - 1) * dt, *popt_line),\
        label='Линейная асимтота', color='m')

ax.set_xlabel(r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$")
ax.set_ylabel(r"$\left\langle r^2\right\rangle \left(\sigma^2\right)$")
ax.set_title("Среднеквадратичное смещение в зависимости от времени")
ax.set_ylim(bottom=-msd[t_range - 2] * 0.07)
ax.legend(loc='best')

fig.set_size_inches(8, 6.3)
fig.savefig("Images/diffusion" + str(counter) + ".png")
fig.savefig("Images/diffusion" + str(counter) + ".svg")
