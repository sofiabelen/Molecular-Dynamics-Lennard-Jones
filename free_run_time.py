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

file_name = 'Data/cross_section'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['x']
y = data['t']

x = np.array(x, dtype=np.float128)
y = np.array(y, dtype=np.float128)

def line(x,a,b):
    return a*x+b

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)

print("mean free run time:",popt[1])

# Write mean free run time to file, to be later read by diff
file_name = 'Data/free_run_time'+str(counter)
f = open(file_name, "a")
f.write(str(popt[1]))
f.close()

plt.style.use('ggplot')
plt.figure(figsize=(20,12))
font_size=25
plt.rc('legend', fontsize=font_size)    # legend fontsize
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

plt.scatter(x,y,label= "Время свободного пробега")
plt.xlabel(r"Время $\left(\sqrt{\frac{\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
plt.ylabel(r"$\tau\left(\sqrt{\frac{\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
plt.title(r"Время свободного пробега", fontsize=font_size)
plt.legend()
plt.savefig("Data/free_run_time"+str(counter)+".png")

plt.legend()
#plt.show()
