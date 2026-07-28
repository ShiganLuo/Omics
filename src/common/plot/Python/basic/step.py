#!/usr/bin/python

###############################################################################
# This script is used for doing the plot of the demographic history of        #
# a random-mating population from a ms command. At the same time, the script  #
# allows to plot (in the same figure) the demographic history infered by the  #
# PSMC software.                                                              #
###############################################################################

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
# Set the values of these global variables
#==============================================================================
# The original ms command:
MS_COMMAND = "./ms 2 100 -t 30000 -r 6000 30000000 -eN 0.01 0.1 -eN 0.06 1 \
                -eN 0.2 0.5 -eN 1 1 -eN 2 2 -p 8 -seeds 1747 45896 23615"

# Path to the output file comming from the PSMC
f1 = "/home/lsg/GWAS/Betta_splendens/PSMC/output/BIM8_combine.psmc"
PSMC_RESULTS = ['BSP2','BIM8','BMA4','BSS4','BSM4']
# Bin size used to generate the imput of PSMC (default is 100)
BIN_SIZE = 100

# Mutation rate per base per generation
MUTATION_RATE = 2.5e-9

# Number of years per generation
GENERAITON_TIME = 1
X_MIN = 1e4

# What plot to do
PLOT_MS = False
PLOT_PSMC_RESULTS = True
#==============================================================================

def ms2fun(ms_command = MS_COMMAND, u = MUTATION_RATE):
    command = ms_command.split(' ')
    N0 = float(command[command.index('-t')+1])/float(command[command.index('-r')+2])/(4*u)
    
    # Getting time and alpha
    size_changes = ms_command.split(' -eN ')
    (t_k, alpha_k) = ([i.split(' ')[0] for i in size_changes[1:]], [j.split(' ')[1] for j in size_changes[1:]])
    
    t0 = min(X_MIN, (GENERAITON_TIME * 4 * N0 * float(t_k[0]))/2)
    # Scalling times and population sizes
    times = [t0] + [GENERAITON_TIME * 4 * N0 * float(i) for i in t_k]
    sizes = [N0] + [N0 * float(i) for i in alpha_k]
    
    times.append(times[-1]*10)
    sizes.append(sizes[-1])
    
    return (times, sizes)

def psmc2fun(filename=f1, s=BIN_SIZE, u=MUTATION_RATE, ith=30):
    #ith represent iteration order,
    a = open(filename, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[ith]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

# if __name__ == "__main__":
    
fig = plt.figure(figsize=(18, 10), dpi=300) 
ax:plt.Axes = fig.subplots() 

if PLOT_MS:
    (real_times, real_sizes) = ms2fun(MS_COMMAND, MUTATION_RATE)
    ax.step(real_times, real_sizes, where='post', linestyle='-', color='k', label = "Real history")

if PLOT_PSMC_RESULTS:
    # for i in range(len(PSMC_RESULTS)):
    #     # PSMC_RESULTS = ["/home/lsg/GWAS/Betta_splendens/PSMC/output/BSP2.psmc","/home/lsg/GWAS/Betta_splendens/PSMC/output/BIM8.psmc","/home/lsg/GWAS/Betta_splendens/PSMC/output/BMA4.psmc"]
    #     # psmc_path = "/home/lsg/GWAS/Betta_splendens/PSMC/output/" + PSMC_RESULTS[i] + ".psmc"
    #     psmc_path = "/home/lsg/GWAS/Betta_splendens/PSMC/output/BSP2_combine.psmc"
    #     color = {'BSP2':'r','BIM8':'#FF7F50','BMA4':'b','BSS4':'#DB7093','BSM4':'#00CED1'}
    #     print(psmc_path)
    #     (estimated_times, estimated_sizes) = psmc2fun(psmc_path, BIN_SIZE, MUTATION_RATE,30) 
    #     ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color=color[PSMC_RESULTS[i]],linewidth=1)
    #     # for j in range(30):
    #     #     (estimated_times, estimated_sizes) = psmc2fun(psmc_path, BIN_SIZE, MUTATION_RATE,j+1)    
    #     #     if i==29:
    #     #         ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color=color[PSMC_RESULTS[i]],linewidth=10)
    #     #     else:
    #     #         ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color=color[PSMC_RESULTS[i]],alpha=0.6,linewidth = 0.2)
    # combine every 1th
    indir = "/home/lsg/GWAS/Betta_splendens/PSMC/output"
    color1 = ['#FF0000','#00ff00','#0000ff','#ff00ff','#ffff00']
    color2 = ['#ff7f7f','#7fff7f','#7f7fff','#ff7fff','#ffff7f']
    for h in range(len(PSMC_RESULTS)):
        psmc_path = f'{indir}/{PSMC_RESULTS[h]}_combine.psmc'
        for i in range(100):
            j = (i+1)*30
            (estimated_times, estimated_sizes) = psmc2fun(psmc_path, BIN_SIZE, MUTATION_RATE,j)
            if j==30:
                ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color=color1[h],linewidth=2)
            elif sum(estimated_sizes)/len(estimated_sizes) == estimated_sizes[0]:
                print(i)
            else:
                ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color=color2[h],alpha=0.5,linewidth =0.5)
legend_elements = [Patch(facecolor=color1[i], label=PSMC_RESULTS[i]) for i in range(len(PSMC_RESULTS))]
ax.legend(handles=legend_elements, loc='upper right', frameon=False)        
    #single 1~30
    # for i in range(30):
    #     psmc_path = "/home/lsg/GWAS/Betta_splendens/PSMC/output/BSP2.psmc"
    #     (estimated_times, estimated_sizes) = psmc2fun(psmc_path, BIN_SIZE, MUTATION_RATE,i+1)    
    #     if i==0:
    #         ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='r',linewidth=0.7)
    #     else:
    #         ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='r',alpha=0.6,linewidth = 0.2)
    # 创建自定义图例项
    # legend_elements = [Patch(facecolor=color[key], label=key) for key in PSMC_RESULTS]
    # 添加图例
    # ax.legend(handles=legend_elements, loc='upper right', frameon=False)
ax.set_xlabel(f"years(μ={MUTATION_RATE},g={GENERAITON_TIME})")
ax.set_ylabel(r'Effective population size($\times10^{4}$)')
def format_func(value, tick_number):
    return f'{value / 1e4:.0f}'

ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_func))
#set scientific number
# ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
# ax.legend(loc='upper right', frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.grid(True)
# ax.set_xlim(X_MIN, X_MAX)
# ax.set_ylim(Y_MIN, Y_MAX)
ax.set_xscale('log')
# ax.set_yscale('log')
fig.savefig("./all.jpg",bbox_inches='tight')
# plt.show()