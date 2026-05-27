import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib as mpl

# Taken from: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
def heatmap(xs, yaxis, xaxis, title, vmin, vmax):
    fig, ax = plt.subplots()
    im = ax.imshow(xs, vmin=vmin, vmax=vmax, cmap="Blues", aspect='0.5')

    ax.set_xlabel('Branching Factor', fontsize='x-large')
    ax.set_ylabel('Number Of Bodies', fontsize='x-large')

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.set_ticks(ticks=[vmin, 1, vmax], labels=[vmin, 1, vmax])
    cbar.ax.set_ylabel(r'Row isolated deviation from BF = 1, ($\frac{bf_i}{bf_1}$)', rotation=-90, va="bottom",fontsize='x-large')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(range(len(xaxis)), labels=xaxis,
                  rotation=45, ha="right", rotation_mode="anchor", fontsize="13")
    ax.set_yticks(range(len(yaxis)), labels=yaxis, fontsize="13")

    # Loop over data dimensions and create text annotations.
    for i in range(len(yaxis)):
        for j in range(len(xaxis)):
            text = ax.text(j, i, xs[i, j],
                           ha="center", va="center", color="#D14F20", fontsize="14" )

    ax.set_title(title, fontsize='xx-large', weight='bold')
    fig.tight_layout()
    plt.show()

def barPlot(data, speedup, xs, title, y2_text, speedup_text ):
    fig, ax = plt.subplots(layout='constrained')
    xticks = ["N = " + str(a) for a in xs]

    x = np.arange(len(xs))  # the label locations
    width = 0.45  # the width of the bars
    multiplier = 0
    for attribute, measurement in data.items():
        offset = width * multiplier
        rects = ax.bar(x + offset , measurement, width, label=attribute)
        ax.bar_label(rects, labels=[str(mea) for mea in measurement],  padding=3) #, fmt=lambda x, p: format(int(x), ','))
        multiplier += 1
    ax.set_xticks(x + width, xs)

    ax.set_xlabel('Number Of Bodies', fontsize='x-large')
    ax.set_ylabel('Runtime in µ seconds', fontsize='x-large')

    ax2 = ax.twinx()
    ax2.plot(x, speedup, color='#1F77B4', marker='s', label=speedup_text)

    ax2.set_ylabel(y2_text, fontsize='x-large')
    ax2.tick_params(axis='y')
    line1, = ax2.plot([0, len(xs)], [1,1], color ="#FF7F0E" , label="Speedup = 1", linestyle='--')
    # fig.legend(loc='outside upper right', fontsize=12)
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left')

    ax.set_title(title, fontsize='xx-large', weight='bold')
    plt.show()

def barPlotForScans(ys_gpu, ys_cpu, ys3, ys4, xs, title, ):
    fig, ax = plt.subplots(layout='constrained')
    data = {
            'Simulated Scan' : ys4,
            'Blocked 64' : ys_gpu,
            'Work Efficient Scan' : ys_cpu,
            'Native Scan' : ys3,
            }
    xticks = ["N = " + str(a) for a in xs]

    x = np.arange(len(xs))  # the label locations
    width = 0.23  # the width of the bars
    multiplier = 0
    for attribute, measurement in data.items():
        offset = width * multiplier
        rects = ax.bar(x + offset , measurement, width, label=attribute)
        ax.bar_label(rects, labels=measurement,  padding=3)
        multiplier += 1
    ax.set_xticks(x + width, xs)

    ax.set_xlabel('Number Of Bodies', fontsize='x-large')
    ax.set_ylabel('Runtime in µ seconds', fontsize='x-large')

    speedup_1 = ys3 / ys_gpu
    speedup_2 = ys3 / ys_cpu
    speedup_3 = ys3 / ys4
    ax2 = ax.twinx()
    ax2.plot(x, speedup_3, color = "#1F77B4", marker ='s', label=r'Speedup: $ \frac{scan}{simulated\ scan}$')
    ax2.plot(x, speedup_1, color = "#FF7F0E", marker ='s', label=r'Speedup: $ \frac{scan}{blocked\ scan\ 64}$')
    ax2.plot(x, speedup_2, color = "#2CA02C", marker ='s', label=r'Speedup: $ \frac{scan}{work\ efficient\ scan}$')

    ax2.set_ylabel('Speedup compared to native scan', fontsize='x-large')
    ax2.tick_params(axis='y')
    line1, = ax2.plot([x[0], x[-1]+1], [1,1], label="Speedup = 1", linestyle='--', color='Red')
    ax2.set(ylim=(0,3) )
    #ax.legend(loc='upper left')
    #ax2.legend(loc='upper left')
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left')
    #fig.legend(loc='upper left')


    ax.set_title(title, fontsize='xx-large', weight='bold')
    plt.show()

def lineGraph(a_gpu, a_cpu, b_gpu, b_cpu, a_xs, b_xs, title, ):
    fig, ax = plt.subplots(layout='constrained')
    x = np.arange(len(a_xs))  # the label locations
    width = 0.23  # the width of the bars
    ax.plot(x, a_gpu, color = "#1F77B4", marker ='s', label=r'Speedup: $ \frac{normal\ RNEA\ gpu}{optimized\ ds\ RNEA\ gpu}$')
    ax.plot(x, a_cpu, color = "#FF7F0E", marker ='s', label=r'Speedup: $ \frac{normal\ RNEA\ cpu}{optimized\ ds\ RNEA\ cpu}$')
    ax.plot([x[0], x[-1]], [1,1], label="Speedup = 1", linestyle='--', color='Red')

    ax.set_xlabel('RNEA: Number Of Bodies', fontsize='x-large')
    ax.set_ylabel('Speedup compared to unoptimized data structure', fontsize='x-large')
    ax.set_ylim(0,6)
    ax.set_xticks(x , a_xs)

    x = np.arange(len(b_xs))  # the label locations
    ax2 = ax.twinx()
    ax2.plot(x, b_gpu, color = "Green", marker ='D', linestyle="dotted",label=r'Speedup: $ \frac{normal\ gpu}{optimized\ ds\ CRBA\ gpu}$')
    ax2.plot(x, b_cpu, color = "#5662F6", marker ='D', linestyle="dotted", label=r'Speedup: $ \frac{normal\ cpu}{optimized\ ds\ CRBA\ cpu}$')
    ax2.set_ylim(0,6)
    ax2.set_yticks([])

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=13)
    ax2 = ax.secondary_xaxis("top", )
    ax2.set_xlabel('CRBA: Number Of Bodies', fontsize='x-large')
    ax2.set_xticks(x , b_xs)

    ax.set_title(title, fontsize='xx-large', weight='bold',y=1.08)
    plt.gcf().set_size_inches(8, 6.8)
    fig.tight_layout()
    plt.savefig("optimized_ds.png", dpi=300)
    plt.show()


def result_to_img(xs, dimentions):
    rnea_img = np.resize(np.array(xs), dimentions)
    rnea_img2 = []
    for x in rnea_img:
        bf1_res = x[0]
        rnea_img2.append([round(a/bf1_res, 4) for a in x])
    rnea_img2 = np.flip(np.array(rnea_img2), 0)
    return rnea_img2

def result_to_img_speedup(xs, dimentions):
    rnea_img = np.resize(np.array(xs), dimentions)
    rnea_img2 = []
    for x in rnea_img:
        bf1_res = x[0]
        rnea_img2.append([round(bf1_res/a, 4) for a in x])
    rnea_img2 = np.flip(np.array(rnea_img2), 0)
    return rnea_img2

rnea_vtree_optimal = np.array([ 657, 656, 657, 660, 668, 668, 669, 668, 751, 765, 767, 763, 1289, 1292, 1289, 1285, 1551, 1639, 1636, 1619, 11318, 14538, 14464, 12730, 21927, 28998, 28610, 24787, 43591, 58851, 57619, 49411 ])
rnea_vtree_optimal_bf1 = rnea_vtree_optimal[0::4]

rnea_vtree_OLD = np.array([ 1902, 1907, 1904, 1903, 1942, 1942, 1943, 1942, 2106, 2112, 2111, 2107, 3316, 3327, 3327, 3320, 7699, 7854, 7839, 7763, 62694, 65913, 65756, 64058, 122196, 129172, 128647, 124869, 242553, 257568, 255818, 247564])
rnea_vtree_old_bf1 = rnea_vtree_OLD[0::4]

rnea_cpu_OLD = np.array([ 4, 4, 4, 3, 29, 28, 28, 28, 281, 277, 278, 277, 2849, 2818, 2818, 2824, 36121, 35801, 35892, 35835, 348002, 345531, 344974, 344668, 688393, 1370050 ])
rnea_cpu_old_bf1 = np.append(rnea_cpu_OLD[0::4], rnea_cpu_OLD[-1])

rnea_cpu_optimal  = np.array([ 3, 3, 5, 5, 19, 18, 18, 19, 165, 167, 165, 165, 1654, 1716, 1664, 1662, 17026, 17070, 16990, 16972, 184429, 181903, 181329, 181503, 397342, 396531, 394820, 393869, 787389, 783728, 784591, 782330])
rnea_cpu_optimal_bf1 = rnea_cpu_optimal[0::4]

rnea_img = result_to_img(rnea_vtree_OLD, (8,4))

#OPTIMAL IS WORSE THAN NON OPTIMAL DS
crba_cpu_optimal =np.array( [9, 9, 12, 8, 7, 370, 369, 176, 131, 74, 33980, 14520, 3391, 2264, 1060, 850975, 154765, 39498, 27364, 17395, 3400138, 419368, 121420, 92538, 63533, 13568473, 1117159, 371170, 297496, 227446, 54322772, 2990806, 1160995, 996526, 838272 ])
crba_cpu_optimal_bf1 = crba_cpu_optimal[0::5]
crba_cpu_optimal_bf1_01 = crba_cpu_optimal[1::5]
crba_cpu_bf2 = crba_cpu_optimal[4::5]

#OLD IS FASTEST!
crba_cpu = np.array([8, 8, 9, 8, 7, 246, 246, 135, 109, 75, 19896, 8736, 2403, 1740, 1017, 511348, 107341, 32332, 23973, 17078, 2008333, 318134, 106074, 84936, 63530, 8079316, 884813, 333416, 282168, 227380, 32462912, 2466171, 1081360, 960038, 849748])
crba_cpu_old_bf1 = crba_cpu[0::5]
crba_cpu_old_bf1_01 = crba_cpu[1::5]
crba_cpu_old_bf1_2 = crba_cpu[3::5]
crba_cpu_old_bf2 = crba_cpu[4::5]

crba_vtree_optimal = np.array([1095, 1086, 1089, 1091, 1092, 1159, 1160, 1161, 1154, 1155, 1250, 1225, 1218, 1278, 1258, 3949, 2456, 2223, 2196, 2175, 13213, 4182, 3564, 3506, 3447, 50620, 7595, 5831, 5708, 5598, 198446, 15398, 10481, 10188, 9903])
crba_gpu_bf1 = crba_vtree_optimal[0::5]
crba_gpu_bf1_01 = crba_vtree_optimal[1::5]
crba_gpu_bf1_2 = crba_vtree_optimal[3::5]
crba_gpu_bf2 = crba_vtree_optimal[4::5]

speedup_crba_bf1 = crba_cpu_old_bf1 / crba_gpu_bf1
speedup_crba_bf1_01 = crba_cpu_old_bf1_01 / crba_gpu_bf1_01
speedup_crba_bf1_2 = crba_cpu_old_bf1_2 / crba_gpu_bf1_2
speedup_crba_bf2 = crba_cpu_old_bf2 / crba_gpu_bf2

data_crba_bf1 = {"vtree gpu" : crba_gpu_bf1, "sequential cpu" : crba_cpu_old_bf1 }
data_crba_bf1_01 = {"vtree gpu" : crba_gpu_bf1_01, "sequential cpu" : crba_cpu_old_bf1_01 }
data_crba_bf1_2 = {"vtree gpu" : crba_gpu_bf1_2, "sequential cpu" : crba_cpu_old_bf1_2 }
data_crba_bf2 = {"vtree gpu" : crba_gpu_bf2, "sequential cpu" : crba_cpu_old_bf2 }

crba_vtree_old = np.array([3593, 3596, 3597, 3632, 3650, 3760, 3758, 3744, 3742, 3744, 4158, 4105, 4057, 4133, 4064, 8568, 6077, 5703, 5676, 5653, 21598, 8377, 7289, 7336, 7296, 73350, 13336, 10897, 10722, 10614, 288754, 25618, 19354, 18991, 18604 ])
crba_vtree_old_bf1 = crba_vtree_old[0::5]
crba_vtree_old_bf1_01 = crba_vtree_old[1::5]

relative_speedup_g_rnea = rnea_vtree_old_bf1 / rnea_vtree_optimal_bf1
relative_speedup_c_rnea = rnea_cpu_old_bf1 / rnea_cpu_optimal_bf1
relative_speedup_g_crba = crba_vtree_old_bf1_01 / crba_gpu_bf1_01
relative_speedup_c_crba = crba_cpu_old_bf1_01 / crba_cpu_optimal_bf1_01

blocked_64     = np.array([609, 698, 1128, 1599, 11046, 21353, ])
work_efficient = np.array([833, 1179, 2009, 4383, 27491, 54540, ])
normal_scan    = np.array([503, 540, 825, 3601, 28978, 54982, ])
simulated_scan = np.array([ 478, 600, 887, 1403, 10040, 19715])

scan_ns = np.array([100, 1000, 10000, 100000, 1000000, 2000000])

blocked_va_512 = np.array([ 206, 235, 362, 1430, 3862, 7224, 10808])
va_unfold = np.array([52, 73, 73, 127, 1617, 3339, ])

crba_img = result_to_img(crba_vtree_optimal, (7,5))

rnea_ns = np.flip(np.array([10, 100, 1000, 10000, 100000, 1000000, 2000000, 4000000]))
rnea_bfs = [1, 1.5 , 2, 1000]
crba_ns = np.flip(np.array([10, 100, 1000, 5000, 10000, 20000, 40000]))
crba_bfs = [1, 1.01, 1.1, 1.2, 2]
 
speedup_rnea = rnea_cpu_optimal_bf1 / rnea_vtree_optimal_bf1
# data_rnea = {"vtree gpu" : rnea_vtree_optimal_bf1, "sequential cpu" : rnea_cpu_optimal_bf1}
speedup_text_rnea = r'Speedup: $\frac{sequential\ cpu}{vtree\ gpu}$'
# barPlot(data_rnea, speedup_rnea, np.flip(rnea_ns), "Performance of RNEA", "Speedup of GPU implementation", speedup_text_rnea)

# heatmap(rnea_img[:7], rnea_ns[:7], [1, 1.5 , 2, 1000], 'RNEA Heatmap', 0.7, 1.3)

# scan_speedup = blocked_64/va_unfold
# speedup_text_scan = r'Speedup: $\frac{Velocity\ Recurrence}{Vector\ Vector\ Addition}$'
# data_scan = {"Mapped scan for vector vector addition" :va_unfold, "Blocked 64 for velcocity recurrence": blocked_64 , }
# barPlot(data_scan, scan_speedup, scan_ns, "Runtime Difference of Scan Operators in rootfix", "Speedup compared to recurrence operator", speedup_text_scan)

# barPlotForScans(blocked_64, work_efficient, normal_scan, simulated_scan, scan_ns, "Performance of Scans")

# heatmap(crba_img, crba_ns, crba_bfs, 'CRBA Heatmap', 0.0, 1)
#
# barPlot(data_crba_bf1, speedup_crba_bf1, np.flip(crba_ns), "Performance of CRBA: Branching Factor = 1", "Speedup of GPU implementation", speedup_text_rnea)
# barPlot(data_crba_bf1_01, speedup_crba_bf1_01, np.flip(crba_ns), "Performance of CRBA: Branching Factor = 1.01", "Speedup of GPU implementation", speedup_text_rnea)
# barPlot(data_crba_bf1_2, speedup_crba_bf1_2, np.flip(crba_ns), "Performance of CRBA: Branching Factor = 1.2", "Speedup of GPU implementation", speedup_text_rnea)
# barPlot(data_crba_bf2, speedup_crba_bf2, np.flip(crba_ns), "Performance of CRBA: Branching Factor = 2", "Speedup of GPU implementation", speedup_text_rnea)

lineGraph(relative_speedup_g_rnea[:7], relative_speedup_c_rnea[:7], relative_speedup_g_crba, relative_speedup_c_crba,  np.flip(rnea_ns[:7]), np.flip(crba_ns),   "Speedup Using Optimized Data Structures")
