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
y = data['x']

popt, pcov = curve_fit(f = line,xdata= x, ydata=y,)

#color1 = "#FF5722"
#color2 = "#536DFE"
#color3 = "#4CAF50"

D = popt[0]/(6.0)

plt.style.use('ggplot')
plt.figure(figsize=(20,10))
plt.title("Среднеквадратичное смещение в зависимость от времени")
plt.xlabel(r"Время $\left(t\right)$")
plt.ylabel(r"Среднеквадратичное смещение $\left(\mathring{A}\right)$")
plt.plot(x, line(x, *popt))
plt.scatter(x,y,label= 'D={0:f}'.format(D))
plt.legend()
plt.savefig("diff.png")
plt.show()
