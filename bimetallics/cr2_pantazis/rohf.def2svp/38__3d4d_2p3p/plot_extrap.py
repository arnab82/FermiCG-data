import os
import sys
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib import colors as mcolors

import scipy

file = os.path.splitext(sys.argv[1])[0]

print(" Working with file: ", file)


np.set_printoptions(suppress=True, precision=6, linewidth=1500)

data_path = file+'.csv'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    headers = next(reader)
    data = list(reader)

conversion = 219474.63 # cm-1
conversion = 1 # hartrees
conversion = 1000 # mH 

energy_var = {}
energy_pt2 = {}
num_thresh = 6
        
n_extrap_points = int((len(data[0])-1)/2)
print(n_extrap_points)

for i in range(len(data)):
    if  i > 0:
        energy_var[i] = np.array([float(a)*conversion for a in data[i][1:n_extrap_points+1]])
        energy_pt2[i] = np.array([float(a)*conversion for a in data[i][n_extrap_points+1:2*n_extrap_points+1]])


emin = 0
for i in energy_pt2:
    emin = min(emin, np.min(energy_pt2[i]))
for i in energy_pt2:
    energy_pt2[i] -= emin
    energy_var[i] -= emin

print(" Energy shift: %12.8f" %emin)

blue = '#3E6D9C'
orange = '#FD841F'
red_orange = '#E14D2A'
dark_blue = '#001253'

cb = ['#000000']
cb.extend([i for i in plt.rcParams['axes.prop_cycle'].by_key()['color']])

fig, ax = plt.subplots()

#this sets the number of decimal points on axis energies
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
extrap = []


#set your x and y axis limits

xmax = 0 
xmin = -10

ymax = 10 
ymin = 0.0
#for s in energy_var:
#    ymax = max(np.max(energy_var[s]), ymax)
#    ymin = min(np.min(energy_pt2[s]), ymin)



for key in energy_var:
    x = energy_pt2[key] - energy_var[key]
    z = energy_pt2[key]
    y = energy_var[key]
    

    m, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    m2, b2, r_value, p_value, std_err = scipy.stats.linregress(x, z)
    
    extrap.append(b)

    plt.rcParams.update({'font.size': 10})

    ymin = min(ymin, b)
    ymax = max(ymax, m*xmin+b)

    ax.plot(x, y, marker='o',linestyle='-' ,markersize=8, color = cb[key])
    ax.plot(x, z, marker='x',linestyle=' ' , markersize=8, color = cb[key])

    x2 = np.array([-1,0])*conversion
    line = m*x2+b
    ax.plot(x2, line, alpha=1.0, color = cb[key], linestyle='-', linewidth=1.5)
    line = m2*x2+b                                     
    ax.plot(x2, line, alpha=0.5, color = cb[key], linestyle='--', linewidth=1.5)
    print("Extrapolated Result: %14.8f"% ((b+emin)/conversion))
    print("R^2                : %14.8f"% r_value)
    print("Var root",key,y)
    print("PT  root",key,z)
    #print("DIFFF",x)
    #print("color",cb[key])


ymin = ymin - 1
print("x: ", (xmin, xmax))
print("y: ", (ymin, ymax))
ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)

ax.set_xlabel('$\Delta$E$_{PT2}$ (mH) ')
ax.set_ylabel('Energy (mH) ')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
#ax.set_yticklabels([])

#ax.legend()
ax.set_title('Absolute Energy, $E(S)$ (mH) ')
#ax.set_title("Tetracence Tetramere \n (40o, 40e) \n ε = {0.0005, 0.0007, 0.001, 0.005, 0.01} \n 31 roots")
#do some legend stuff or comment out
#black_patch = mpatches.Patch(color=dark_blue, label='Ground')
#blue_patch = mpatches.Patch(color=blue, label='Triplets')
#green_patch = mpatches.Patch(color=orange, label='Singlets')
#red_patch = mpatches.Patch(color=red_orange, label='Biexcitons')
#ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='center right')
#ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='upper center')
#ax.legend()

fig = plt.gcf()
fig.set_size_inches(4.5,4.5)
fig.savefig(file+'_extrap.pdf', dpi=300, bbox_inches='tight')

#plt.show()

print(extrap)
