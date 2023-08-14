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
import matplotlib.lines as mlines


file = os.path.splitext(sys.argv[1])[0]
print(" Working with file: ", file)

nroots = 4
np.set_printoptions(suppress=True, precision=6, linewidth=1500)

data_path = file+'.csv'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    headers = next(reader)
    data = list(reader)

conversion = 1000 # mH 
conversion = 1 # hartrees
conversion = 219474.63 # cm-1

energy_var = {}
energy_pt2 = {}
num_thresh = 6

n_extrap_points = int((len(data[0])-1)/2)
print(n_extrap_points)

for i in range(len(data)):
    if  i > 0:
        energy_var[i] = np.array([float(a)*conversion for a in data[i][1:n_extrap_points+1]])
        energy_pt2[i] = np.array([float(a)*conversion for a in data[i][n_extrap_points+1:2*n_extrap_points+1]])


gaps_var = {} 
gaps_pt2 = {}


for i in range(2,nroots+1):
    gaps_var[i] = energy_var[i] - energy_var[i-1]
    gaps_pt2[i] = energy_pt2[i] - energy_pt2[i-1]




cb = ['#000000']
cb.extend([i for i in plt.rcParams['axes.prop_cycle'].by_key()['color']])

fig, ax = plt.subplots()

#this sets the number of decimal points on axis energies
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
extrap = []


#set your x and y axis limits
ymax = -1e10 
ymin = 0.0
for s in gaps_var:
    ymax = max(np.max(gaps_var[s]), ymax)
    ymin = min(np.min(gaps_pt2[s]), ymin)



ymax = 0 
ymin = 250 
xmax = 0 
xmin = 30.0

for key in gaps_var:
    print(key, " ", gaps_var[key], gaps_pt2[key])
    x = gaps_pt2[key] - gaps_var[key]
    z = gaps_pt2[key]
    y = gaps_var[key]


   
    ax.plot(x, y, marker='o',linestyle='-' ,markersize=10, color = cb[key])
    ax.plot(x, z, marker='x',linestyle=' ' , markersize=10, color = cb[key])
    
    continue
    m, b = np.polyfit(x, y, 1)
    m2,b = np.polyfit(x, z, 1)
    extrap.append(b)

    plt.rcParams.update({'font.size': 10})

    print("b: ", b)
    ymin = min(ymin, b)
    ymax = max(ymax, m*xmin+b)


    x2 = np.array([-1,0])*0.9
    line = m*x2+b
    ax.plot(x2, line, 'r', alpha=1.0, color = cb[key], linestyle='-', linewidth=1.5)
    line = m2*x2+b                                     
    ax.plot(x2, line, 'r', alpha=0.5, color = cb[key], linestyle='--', linewidth=1.5)
    print("Extrap  %14.8f"%b)
    print("Var root",key,y)
    print("PT  root",key,z)
    print("DIFFF",x)
    print("color",cb[key])

var_marker = mlines.Line2D([], [], color='grey', marker='o', linestyle='None',
                          markersize=10, label='Variational')
pt2_marker = mlines.Line2D([], [], color='grey', marker='x', linestyle='None',
                          markersize=10, label='PT2')


print("x: ", (xmin, xmax))
print("y: ", (ymin, ymax))
ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)

ax.set_xlabel('$\Delta\Delta$E$_{PT2}$ (cm$^{-1}$) ')
ax.set_ylabel('Gap (cm$^{-1}$) ')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
#ax.set_yticklabels([])

ax.legend(handles=[var_marker, pt2_marker], loc='upper left')

#ax.set_title("Tetracence Tetramere \n (40o, 40e) \n ε = {0.0005, 0.0007, 0.001, 0.005, 0.01} \n 31 roots")
ax.set_title('Energy Gaps, $E(S)-E(S-1)$ (cm$^{-1}$) ')
#do some legend stuff or comment out
#black_patch = mpatches.Patch(marker='x', label='PT2')
#blue_patch = mpatches.Patch(color=blue, label='Triplets')
#green_patch = mpatches.Patch(color=orange, label='Singlets')
#red_patch = mpatches.Patch(color=red_orange, label='Biexcitons')
#ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='center right')
#ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='upper center')
#ax.legend()

fig = plt.gcf()
fig.set_size_inches(4.5,4.5)
fig.savefig(file+'_gaps.pdf', dpi=300, bbox_inches='tight')

#plt.show()

print(extrap)
