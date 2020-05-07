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
counter = 22

file_name = 'Data/pv'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['x']
y = data['y']

x = np.array(x, dtype=np.float128)
y = np.array(y, dtype=np.float128)

def line(x,a,b):
    return a*x+b

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)
sigma = np.sqrt(np.diag(pcov))
print("sigma=",sigma[0])

plt.style.use('ggplot')
plt.figure(figsize=(20,10))

font_size=25
plt.figure(figsize=(20,12))
plt.rc('legend', fontsize=font_size)    # legend fontsize
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

plt.scatter(x,y,label= "PV")

#plt.plot(x, line(x, *popt),c='black')
plt.xlabel(r"Время $\left(\sqrt{\frac{\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
plt.ylabel(r"PV", fontsize=font_size)
plt.title(r"Зависимость $PV$ от времени", fontsize=font_size)
plt.legend()
plt.savefig("Data/pv"+str(counter)+".png")
#plt.show()
