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
#counter = 16

file_name = 'Data/temperature'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['x']
y = data['y']

x = np.array(x, dtype=np.float128)
y = np.array(y, dtype=np.float128)

def line(x,a,b):
    return a*x+b

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)
sigma = np.sqrt(np.diag(pcov))
print("sigma temp=",sigma[0])

plt.style.use('ggplot')
plt.figure(figsize=(20,12))
font_size=25
plt.rc('legend', fontsize=font_size)    # legend fontsize
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

#plt.rc('font', size=font_size)          # controls default text sizes
#plt.rc('axes', titlesize=font_size)     # fontsize of the axes title
#plt.rc('axes', labelsize=font_size)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=font_size)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=font_size)    # fontsize of the tick labels
#plt.rc('legend', fontsize=font_size)    # legend fontsize
#plt.rc('figure', titlesize=font_size)

plt.scatter(x,y)
plt.xlabel(r"Время $\left(\sqrt{\frac{\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
plt.ylabel(r"Температура $\left(\frac{\varepsilon}{k_{Б} }\right)$", fontsize=font_size)
plt.title(r"Зависимость температуры от времени", fontsize=font_size)
plt.legend()
plt.savefig("Data/temperature"+str(counter)+".png")
plt.show()
