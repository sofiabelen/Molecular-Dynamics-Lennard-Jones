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
file_name = 'Data/velocity'+str(counter)

with open(file_name) as f:
    tmp = [float(x) for x in next(f).split()] # read first line
    n = int(tmp[0])
    m = int(tmp[1])
    m = m*3
    bins = int(tmp[2])
    #T = tmp[3]
    i=0
    tmp = np.zeros((m,n))
    for line in f: # read rest of lines
        tmp[i]=([float(x) for x in line.split()])
        i = i+1
        if(i == m):
            break

vel = np.zeros(m*n)
max_v=0.0
for i in range(m):
    for j in range(n):
        vel[i*n+j]=tmp[i][j]**2
        max_v=max(max_v,vel[i*n+j])

x = np.linspace(0.0,int(max_v)+1,num = int(bins))
#y = np.zeros(bins)
#for i in range(bins):
#    y[i] = mt.exp(-x[i]/2*T)
#plt.plot(x,y,label = r'$e^{\frac{mv_x^2}{2kT}}$')

plt.style.use('ggplot')

# -------------------- #
p, x = np.histogram(vel, bins=bins)
x = x[:-1] + (x[1] - x[0])/2     # convert bin edges to centers
f = UnivariateSpline(x, p, s=100)
plt.plot(x, f(x))
# -------------------- #

plt.title(r'Распределение скоростей по проекциям')
plt.hist(vel,bins=bins)
plt.yscale('log', nonposy='clip')
plt.xlabel(r"$v_x^2$")
plt.ylabel(r"$\log_{10}(n_v)$")
plt.legend()

plt.savefig("media/hist.png")
plt.savefig("media/hist.svg")

plt.show()
