import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import math as mt
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline

trk = open('Data/counter','r')
counter = int(trk.read()) - 1
#counter = 19

file_name = 'Data/energy'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['x']
y = data['e']
y2 = data['k']
y3 = data['p']

#x = np.array(x, dtype=np.float128)
#y = np.array(y, dtype=np.float128)
#
#def line(x,a,b):
#    return a*x+b
#
#popt, pcov = curve_fit(f = line, xdata = x, ydata = y)
#sigma = np.sqrt(np.diag(pcov))
#print("sigma energy=",sigma[0])
#popt, pcov = curve_fit(f = line, xdata = x, ydata = y2)
#sigma = np.sqrt(np.diag(pcov))
#print("sigma kinetic=",sigma[0])
#popt, pcov = curve_fit(f = line, xdata = x, ydata = y3)
#sigma = np.sqrt(np.diag(pcov))
#print("sigma potential=",sigma[0])

plt.style.use('ggplot')
font_size=35
plt.figure(figsize=(25,14))
plt.rc('legend', fontsize=font_size)    # legend fontsize
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

#plt.scatter(x,y,label=r"$\langle E\rangle$")
#plt.scatter(x,y2,label=r"$\langle K\rangle$")
#plt.scatter(x,y3,label=r"$\langle П \rangle$")

plt.scatter(x,y,label= "Средняя полная энергия")
#plt.plot(x, line(x, *popt),c='black')
plt.xlabel(r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
plt.ylabel(r"Энергия $\left(\varepsilon\right)$", fontsize=font_size)
plt.title(r"Зависимость энергии от времени", fontsize=font_size)
plt.legend()

#plt.show()
plt.savefig("Data/energy"+str(counter)+".png")

plt.scatter(x,y2,label="Средняя кинетическая энергия")
plt.scatter(x,y3,label="Средняя потенциальная энергия")
plt.legend()

plt.savefig("Data/energy_kp"+str(counter)+".png")
plt.savefig("../Lab2/selfD/energy_kp"+str(counter)+".png")
#plt.show()
