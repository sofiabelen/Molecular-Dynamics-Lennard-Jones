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
counter = 54

file_name = 'Data/diff'+str(counter)
data = pd.read_table(file_name,sep = '\s+')
x = data['t']
y = data['x']

x = np.array(x, dtype=np.float128)
y = np.array(y, dtype=np.float128)

popt, pcov = curve_fit(f = line, xdata = x, ydata = y)

#Read mean free run time
file_name = 'Data/free_run_time'+str(counter)
f = open(file_name, "r")
T = float((f.read()))

file_name = 'Data/diffusion_mol'+str(counter)
f = open(file_name, "r")
D_mol = float((f.read()))

plt.style.use("ggplot")

font_size=30
plt.figure(figsize=(17,15))
plt.rc('legend', fontsize=30)    # legend fontsize
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

plt.scatter(x,y, label=r"$\left\langle r^2\right\rangle(t)$")

slinear = 0
for i in range(len(x)):
    if(x[i] >= 4.0):
        slinear = i
        break

eparab = 0
for i in range(len(x)):
    if(x[i] >= 1.0):
        eparab = i
        break
    
x2 = np.delete(x, np.s_[eparab:int(len(x))])
y2 = np.delete(y, np.s_[eparab:int(len(y))])

x3 = np.delete(x, np.s_[:slinear])
y3 = np.delete(y, np.s_[:slinear])

#x = np.delete(x, np.s_[:int(len(x)/4)])
#y = np.delete(y, np.s_[:int(len(y)/4)])

popt, pcov = curve_fit(f = line, xdata = x3, ydata = y3)
popt2, pcov2 = curve_fit(f = parabola, xdata = x2, ydata = y2)

#sigma = np.sqrt(np.diag(pcov))
#print("sigma=",sigma)
#print(pcov)

D = popt[0]/2.0
print("slope: D=", D);
print("free run time: D=", D_mol);

#plt.xlim(left=1.6)
#plt.xlim(right=1.5)
#plt.ylim(bottom=-0.5)
#plt.ylim(top=1.0)
#plt.xscale("log")
#plt.yscale("log")

#plt.loglog(x,y)
#plt.xlim(left=0.1)
#plt.ylim(bottom=0.01)

#plt.xlim(right=5)
#plt.ylim(bottom=-0.5)
#plt.ylim(top=8.0)

plt.ylim(top=9.8)

x1 = np.log(5)/np.log(10)
x2 = np.log(4)/np.log(10)
y1 = np.log(line(5, *popt))/np.log(10)
y2 = np.log(line(4, *popt))/np.log(10)

x3 = np.log(0.1)/np.log(10)
x4 = np.log(1)/np.log(10)
y3 = np.log(parabola(0.1, *popt2))/np.log(10)
y4 = np.log(parabola(1, *popt2))/np.log(10)

m1 = (y2-y1)/(x2-x1)
b1 = - x1*m1 + y1

m2 = (y4-y3)/(x4-x3)
b2 = - x3*m2 + y3

xi = (b2-b1)/(m1-m2)
yi = m1*xi + b1

print("free run time (extrapol) = ",10**xi)

def ex_line(x, a, b):
   return 10**(a * np.log(x)/np.log(10) + b)

#plt.plot(10**xi, 10**yi, color='black')
#plt.plot(x, ex_line(x, m1, b1),label= 'Экстраполяция линейного участка',color="#2CA02C")
#plt.plot(x, ex_line(x, m1, b1),color="#2CA02C")
#plt.plot(x, ex_line(x, m2, b2),label= 'Экстраполяция параболического участка',color="#1F77B4")
#plt.plot(10**xi,10**yi,'ro', color='black')

plt.plot(x, line(x, *popt),label= 'Линейная асимтота',color="#2CA02C")
plt.plot(x, parabola(x, *popt2),label= 'Параболическая асимтота',color="#1F77B4")

#plt.plot(x, line(x, *popt),label= 'D={0:f}'.format(D),color="black")

#plt.axvline(x=T, color='k', label=r"Время свободного пробега $\tau=$"+"%.2f"%T)

plt.xlabel(r"Время $\left(\sqrt{\frac{m\sigma^2}{\varepsilon} }\right)$", fontsize=font_size)
plt.ylabel(r"$\left\langle r^2\right\rangle \left(\sigma^2\right)$", fontsize=font_size)
plt.title("Среднеквадратичное смещение в зависимость от времени", fontsize=font_size+8)
plt.legend()
plt.savefig("Data/diff"+str(counter)+".png")
plt.savefig("../Lab2/selfD/diff"+str(counter)+".png")
#plt.show()
