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

file_name = 'Data/diff'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['t']
y = data['x']

x = np.array(x, dtype=np.float128)
y = np.array(y, dtype=np.float128)

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)

#color1 = "#FF5722"
#color2 = "#536DFE"
#color3 = "#4CAF50"



plt.style.use('ggplot')
plt.figure(figsize=(20,10))

plt.scatter(x,y)

popt, pcov = curve_fit(f = parabola, xdata = x, ydata = y)
#plt.plot(x, parabola(x, *popt))

print(len(x))
x = np.delete(x, np.s_[:int(len(x)/2)])
y = np.delete(y, np.s_[:int(len(y)/2)])
print(len(x))

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)
D = popt[0]/(6.0)
plt.plot(x, line(x, *popt),label= 'D={0:f}'.format(D))

plt.xlabel(r"Время $\left(t\right)$")
plt.ylabel(r"Среднеквадратичное смещение $\left(\sigma^2\right)$")
plt.title("Среднеквадратичное смещение в зависимость от времени")

plt.legend()
plt.savefig("Data/diff"+str(counter)+".png")
plt.show()
