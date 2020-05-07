import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import math as mt
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline

def line(x,a,b):
    return a*x+b

def parabola(x,a,b,c):
    return a*x*x+b*x+c

trk = open('Data/counter','r')
counter = int(trk.read()) - 1
counter = 17

file_name = 'Data/diff'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['t']
y = data['x']

x = np.array(x, dtype=np.float128)
y = np.array(y, dtype=np.float128)

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)

plt.style.use("ggplot")

font_size=25
plt.figure(figsize=(20,12))
plt.rc('legend', fontsize=20)    # legend fontsize
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

plt.scatter(x,y)

x2 = np.delete(x, np.s_[int(len(x)/7):int(len(x))])
y2 = np.delete(y, np.s_[int(len(y)/7):int(len(y))])

x3 = np.delete(x, np.s_[:int(len(x)/4)])
y3 = np.delete(y, np.s_[:int(len(y)/4)])

#x = np.delete(x, np.s_[:int(len(x)/4)])
#y = np.delete(y, np.s_[:int(len(y)/4)])

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)

sigma = np.sqrt(np.diag(pcov))
print("sigma=",sigma)
print(pcov)

D = popt[0]*0.5

plt.xlim(left=1.6)
plt.xlim(right=10)
plt.ylim(bottom=0.8)
plt.xscale("log")
plt.yscale("log")

popt, pcov = curve_fit(f = line, xdata = x3, ydata = y3)
plt.plot(x, line(x, *popt),label= 'Линейное фитирование',color="black")

popt, pcov = curve_fit(f = parabola, xdata = x2, ydata = y2)
plt.plot(x, parabola(x, *popt),label= 'Параболическое фитирование',color="grey")

#plt.plot(x, line(x, *popt),label= 'D={0:f}'.format(D),color="black")
plt.xlabel(r"Время $\left(\sqrt{\frac{\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
plt.ylabel(r"Среднеквадратичное смещение $\left(\sigma^2\right)$", fontsize=font_size)
plt.title("Среднеквадратичное смещение в зависимость от времени", fontsize=font_size)
plt.legend()
plt.savefig("Data/diff_log"+str(counter)+".png")
plt.show()
