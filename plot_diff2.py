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

data = pd.read_table('diff',sep = '\s+')
x = data['t']
y1 = data['x']
y2 = data['y']
y3 = data['z']

poptx, pcovx = curve_fit(f = line,xdata= x, ydata=y1,)
popty, pcovy = curve_fit(f = line,xdata= x, ydata=y2,)
poptz, pcovz = curve_fit(f = line,xdata= x, ydata=y3,)

color1 = "#FF5722"
color2 = "#536DFE"
color3 = "#4CAF50"

plt.style.use('ggplot')
plt.title("Среднеквадратичное смещение в зависимость от времени")
plt.xlabel(r"Время $\left(t\right)$")
plt.ylabel(r"Среднеквадратичное смещение $\left(\mathring{A}\right)$")
plt.plot(x, line(x, *poptx),c=color1)
plt.plot(x, line(x, *popty),c=color2)
plt.plot(x, line(x, *poptz),c=color3)
plt.scatter(x,y1,label=r'СКС - x',c=color1)
plt.scatter(x,y2,label=r'СКС - y',c=color2)
plt.scatter(x,y3,label=r'СКС - z',c=color3)
plt.legend()
plt.savefig("diff.png")
plt.show()
