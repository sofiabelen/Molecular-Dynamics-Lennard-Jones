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
dt = param[5]

data = data.reshape((n_steps, n_part, dim))

t_range = int((n_steps-1) / 2)

msd = np.zeros(t_range - 1)

for deltaT in range(1, t_range):
    for d in range(0, dim):
        for i in range(0, n_part):
            for t0 in range(0, t_range):
                msd[deltaT - 1] +=\
                        (data[t0 + deltaT][i][d] - data[t0][i][d])**2
    msd[deltaT - 1] /= float(n_part * t_range * dim)


def line(x, a, b):
    return a*x + b

def parabola(x, a, b, c):
    return a*x**2 + b*x + c

parabola_end = int(t_range * 0.1)
line_start = int(t_range * 0.6)
time = np.arange(0, t_range - 1) * dt

popt_parabola, pcov_parabola = curve_fit(f=parabola,\
        xdata=np.arange(0, parabola_end) * dt,\
        ydata=msd[0:parabola_end])

popt_line, pcov_line = curve_fit(f=line,\
        xdata=np.arange(line_start, t_range - 1) * dt,\
        ydata=msd[line_start:t_range - 1])


sns.set(context='notebook', style='darkgrid')
pal = sns.hls_palette(10, l=.7, h=.5, s=.4)
pal_dark = sns.hls_palette(10, l=.3)
color_scatter = pal[0]
color_line = pal_dark[7]
color_parab = pal_dark[3]

width_scatter = 4.0
width_lines = 3.5
width_point = 5.0

fig, ax = plt.subplots()

ax.scatter(time, msd, color=color_scatter, linewidth=width_scatter)

ax.plot(time, line(time, *popt_line),\
        label='Линейная асимтота', color=color_line,\
        linewidth=width_lines)

ax.plot(time, parabola(time, *popt_parabola),\
        label='Параболическая асимтота', color=color_parab,\
        linewidth=width_lines)

ax.set_xlabel(r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$")
ax.set_ylabel(r"$\left\langle r^2\right\rangle \left(\sigma^2\right)$")
ax.set_title("Среднеквадратичное смещение в зависимости от времени")
ax.set_ylim(bottom=-msd[t_range - 2] * 0.07)
ax.legend(loc='best')

fig.set_size_inches(8, 6.3)
fig.savefig("Images/diffusion" + str(counter) + ".png")
fig.savefig("Images/diffusion" + str(counter) + ".svg")

# --- Log-Log Plot --- #

fig_log, ax_log = plt.subplots()

ax_log.scatter(time, msd, color=color_scatter, linewidth=width_scatter)

ax_log.set_xscale('log')
ax_log.set_yscale('log')

def line_in_log(x, a, b):
   return 10**(a * np.log(x) / np.log(10) + b)

value_line1 = (t_range - t_range * 0.1 - 1.0) * dt
value_line2 = (t_range - t_range * 0.3) * dt
value_parab1 = t_range * 0.01 * dt
value_parab2 = t_range * 0.1 * dt

x1 = np.log(value_line1) / np.log(10)
x2 = np.log(value_line2) / np.log(10)
y1 = np.log(line(value_line1, *popt_line)) / np.log(10)
y2 = np.log(line(value_line2, *popt_line)) / np.log(10)

x3 = np.log(value_parab1) / np.log(10)
x4 = np.log(value_parab2) / np.log(10)
y3 = np.log(parabola(value_parab1, *popt_parabola)) / np.log(10)
y4 = np.log(parabola(value_parab2, *popt_parabola)) / np.log(10)

m1 = (y2 - y1) / (x2 - x1)
b1 = -x1*m1 + y1

m2 = (y4 - y3) / (x4 - x3)
b2 = -x3*m2 + y3

xi = (b2 - b1) / (m1 - m2)
yi = m1*xi + b1

time = np.arange(1, t_range - 1) * dt

ax_log.plot(time, line_in_log(time, m1, b1),\
        label='Экстраполяция линейного участка', color=color_line,\
        linewidth=width_lines)

ax_log.plot(time, line_in_log(time, m2, b2),\
        label='Экстраполяция параболического участка',\
        color=color_parab, linewidth=width_lines)

ax_log.plot(10**xi,10**yi, 'ro',\
        label=r'Время свободного пробега $\tau=$'+'%.2f'%10**xi,\
        color='black', linewidth=width_point)

ax_log.set_xlabel(\
        r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$")
ax_log.set_ylabel(\
        r"$\left\langle r^2\right\rangle \left(\sigma^2\right)$")
ax_log.set_title(\
        "Среднеквадратичное смещение в зависимости от времени")
ax_log.legend(loc='best')

fig_log.set_size_inches(8, 6.3)
fig_log.savefig("Images/diffusion_log" + str(counter) + ".png")
fig_log.savefig("Images/diffusion_log" + str(counter) + ".svg")
