import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import math as mt
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline
import statistics as stats

trk = open('Data/counter','r')
counter = int(trk.read()) - 1
#counter = 37

file_name = 'Data/cross_section'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['x']
t = data['t']
d = data['d']

x = np.array(x, dtype=np.float128)
t = np.array(t, dtype=np.float128)
d = np.array(d, dtype=np.float128)

# Write mean free run time to file, to be later read by diff
file_name = 'Data/free_run_time'+str(counter)
f = open(file_name, "w")
f.write(str(stats.mean(t)))
f.close()

file_name = 'Data/diffusion_mol'+str(counter)
f = open(file_name, "w")
f.write(str(stats.mean(d)))
f.close()

#plt.style.use('ggplot')
#plt.figure(figsize=(20,12))
#font_size=25
#plt.rc('legend', fontsize=font_size)    # legend fontsize
#plt.xticks(fontsize=font_size)
#plt.yticks(fontsize=font_size)
#
#plt.scatter(x,y,label= "Время свободного пробега")
#plt.xlabel(r"Время $\left(\sqrt{\frac{\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
#plt.ylabel(r"$\tau\left(\sqrt{\frac{\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
#plt.title(r"Время свободного пробега", fontsize=font_size)
#plt.legend()
#plt.savefig("Data/free_run_time"+str(counter)+".png")
#
#plt.legend()
#plt.show()
