import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib as mpl

# Taken from: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
def heatmap(xs, yaxis, xaxis, title, vmin, vmax):
    fig, ax = plt.subplots()
    im = ax.imshow(xs, vmin=vmin, vmax=vmax, cmap="Blues")

    ax.set_xlabel('Branching Factor', fontsize='x-large')
    ax.set_ylabel('Number Of Bodies', fontsize='x-large')

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.set_ticks(ticks=[vmin, 1, vmax], labels=[vmin, 1, vmax])
    cbar.ax.set_ylabel(r'Row isolated deviation from BF = 1, ($\frac{bf_i}{bf_1}$)', rotation=-90, va="bottom",fontsize='large')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(range(len(xaxis)), labels=xaxis,
                  rotation=45, ha="right", rotation_mode="anchor", )
    ax.set_yticks(range(len(yaxis)), labels=yaxis)

    # Loop over data dimensions and create text annotations.
    for i in range(len(yaxis)):
        for j in range(len(xaxis)):
            text = ax.text(j, i, xs[i, j],
                           ha="center", va="center", color="Red", fontsize="small" )

    ax.set_title(title, fontsize='xx-large')
    fig.tight_layout()
    plt.savefig(title + ".png", dpi=300)
    plt.show()

def barPlot(ys_gpu, ys_cpu, xs, title, ):
    fig, ax = plt.subplots(layout='constrained')
    data = {
            'Blocked 64 for velocity recurrence' : ys_cpu,
            'Blocked 512 for vector addition' : ys_gpu,
            }
    xticks = ["N = " + str(a) for a in xs]

    x = np.arange(len(xs))  # the label locations
    width = 0.35  # the width of the bars
    multiplier = 0
    for attribute, measurement in data.items():
        offset = width * multiplier
        rects = ax.bar(x + offset , measurement, width, label=attribute)
        ax.bar_label(rects, labels=measurement,  padding=3)
        multiplier += 1
    ax.set_xticks(x + width, xs)

    ax.set_xlabel('Number Of Bodies', fontsize='x-large')
    ax.set_ylabel('Runtime in µ seconds', fontsize='x-large')

    speedup_1 = ys_cpu / ys_gpu
    ax2 = ax.twinx()
    ax2.plot(x, speedup_1, color='#eb7e1fff', marker='s', label=r'Speedup: $\frac{Velocity\_Recurrence}{Vector\_Vector\_Addition}$')

    ax2.set_ylabel('Speedup compared to recurrence operator', fontsize='x-large')
    ax2.tick_params(axis='y', labelcolor='b')
    line1, = ax2.plot([0, len(xs)], [1,1], label="Speedup = 1", linestyle='--')
    fig.legend(loc='outside upper right', fontsize=12)

    ax.set_title(title, fontsize='xx-large', weight='bold')
    plt.show()

def barPlotForScans(ys_gpu, ys_cpu, ys3, xs, title, ):
    fig, ax = plt.subplots(layout='constrained')
    data = {
            'Blocked_64' : ys_gpu,
            'Work Efficient Scan' : ys_cpu,
            'Native Scan' : ys3
            }
    # linegpu, = ax.plot(xs, ys_gpu, marker="D")
    # linecpu, = ax.plot(xs, ys_cpu, marker="D")
    # linegpu.set_label('vtree gpu')
    # linecpu.set_label('sequential cpu')
    xticks = ["N = " + str(a) for a in xs]

    x = np.arange(len(xs))  # the label locations
    width = 0.25  # the width of the bars
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
    ax2 = ax.twinx()
    ax2.plot(x, speedup_1, 'b-s', label=r'Speedup: $ \frac{scan}{blocked\_scan\_64}$')
    ax2.plot(x, speedup_2, 'r-s', label=r'Speedup: $ \frac{scan}{work\_efficient\_scan}$')

    ax2.set_ylabel('Speedup compared to native scan', fontsize='x-large')
    ax2.tick_params(axis='y', labelcolor='b')
    line1, = ax2.plot([0, len(xs)], [1,1], label="Speedup = 1", linestyle='--')
    fig.legend(loc='outside upper right')


    ax.set_title(title, fontsize='xx-large', weight='bold')
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
rnea_vtree_bf1 = rnea_vtree_OLD[0::4]

rnea_cpu_OLD = np.array([ 4, 4, 4, 3, 29, 28, 28, 28, 281, 277, 278, 277, 2849, 2818, 2818, 2824, 36121, 35801, 35892, 35835, 348002, 345531, 344974, 344668, 688393, 1370050 ])

rnea_cpu_optimal  = np.array([ 3, 3, 5, 5, 19, 18, 18, 19, 165, 167, 165, 165, 1654, 1716, 1664, 1662, 17026, 17070, 16990, 16972, 184429, 181903, 181329, 181503, 397342, 396531, 394820, 393869, 787389, 783728, 784591, 782330])
rnea_cpu_optimal_bf1 = rnea_cpu_optimal[0::4]

rnea_img = result_to_img(rnea_vtree, (8,4))

#OPTIMAL IS WORSE THAN NON OPTIMAL DS
crba_cpu_optimal =np.array( [9, 9, 12, 8, 7, 370, 369, 176, 131, 74, 33980, 14520, 3391, 2264, 1060, 850975, 154765, 39498, 27364, 17395, 3400138, 419368, 121420, 92538, 63533, 13568473, 1117159, 371170, 297496, 227446, 54322772, 2990806, 1160995, 996526, 838272 ])
crba_cpu_bf1 = np.array([8, 8, 9, 8, 7, 246, 246, 135, 109, 75, 19896, 8736, 2403, 1740, 1017, 511348, 107341, 32332, 23973, 17078, 2008333, 318134, 106074, 84936, 63530, 8079316, 884813, 333416, 282168, 227380, 32462912, 2466171, 1081360, 960038, 849748])
crba_cpu_bf1 = crba_cpu[0::4]
crba_cpu_bf1 = crba_cpu[1::4]
crba_cpu_bf1 = crba_cpu[1::4]
crba_cpu_bf1 = crba_cpu[2::4]

crba_vtree_optimal = np.array([1095, 1086, 1089, 1091, 1092, 1159, 1160, 1161, 1154, 1155, 1250, 1225, 1218, 1278, 1258, 3949, 2456, 2223, 2196, 2175, 13213, 4182, 3564, 3506, 3447, 50620, 7595, 5831, 5708, 5598, 198446, 15398, 10481, 10188, 9903])

crba_vtree_old = np.array([3593, 3596, 3597, 3632, 3650, 3760, 3758, 3744, 3742, 3744, 4158, 4105, 4057, 4133, 4064, 8568, 6077, 5703, 5676, 5653, 21598, 8377, 7289, 7336, 7296, 73350, 13336, 10897, 10722, 10614, 288754, 25618, 19354, 18991, 18604 ])
crba_gpu_bf1 = crba_vtree[0::4]
crba_gpu_bf1 = crba_vtree[0::4]
crba_gpu_bf1 = crba_vtree[0::4]
crba_gpu_bf1 = crba_vtree[4::4]


blocked_64 = np.array([1114, 1197, 2416, 3524, 21629, 41474, 62510])
work_efficient = np.array([ 625, 933, 2038, 8402, 64390, 131305, 260649])
normal_scan = np.array([937, 1163, 2393, 14639, 129984, 246428, 363029])

blocked_va_512 = np.array([ 206, 235, 362, 1430, 3862, 7224, 10808])

crba_img = result_to_img(crba_vtree, (7,5))

rnea_ns = np.flip(np.array([10, 100, 1000, 10000, 100000, 1000000, 2000000, 4000000]))
rnea_bfs = [1, 1.5 , 2, 1000]
crba_ns = np.flip(np.array([10, 100, 1000, 5000, 10000, 20000, 40000]))
crba_bfs = [1, 1.01, 1.1, 1.2, 2]
 
# print(rnea_vtree_bf1, rnea_cpu_bf1)
# barPlot(rnea_vtree_bf1, rnea_cpu_bf1, np.flip(rnea_ns), "Performance of RNEA")
barPlot(blocked_va_512, blocked_64, np.flip(rnea_ns)[1:], "Runtime Difference of Scan Operators in rootfix")
# barPlotForScans(blocked_64, work_efficient, normal_scan, np.flip(rnea_ns)[1:], "Performance of Scans")
# heatmap(rnea_img[:7], rnea_ns[:7], [1, 1.5 , 2, 1000], 'RNEA Heatmap', 0.96, 1.04)
# heatmap(crba_img, crba_ns, crba_bfs, 'CRBA Heatmap', 0.0, 1)
