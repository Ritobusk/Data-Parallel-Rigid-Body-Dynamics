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
            'sequential cpu' : ys_cpu,
            'vtree gpu' : ys_gpu,
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
    ax2.plot(x, speedup_1, color='#eb7e1fff', marker='s', label=r'Speedup: $\frac{sequential\_cpu}{vtree\_gpu}$')

    ax2.set_ylabel('Speedup of GPU implementation', fontsize='x-large')
    ax2.tick_params(axis='y', labelcolor='b')
    line1, = ax2.plot([0, len(xs)], [1,1], label="Speedup = 1", linestyle='--')
    fig.legend(loc='outside upper right', fontsize=12)

    ax.set_title(title, fontsize='xx-large', weight='bold')
    plt.savefig(title + ".png", dpi=300)
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

rnea_vtree = np.array([1541, 1539, 1539, 1539, 1606, 1606, 1611, 1606, 1794, 1796, 1794, 1799, 3378, 3384, 3386, 3389, 8165, 8265, 8307, 8237, 46632, 48797, 48648, 47295, 89426, 94340, 93850, 90875, 133928, 141491, 140606, 136116 ])
rnea_vtree_bf1 = rnea_vtree[0::4]

rnea_cpu  = np.array([5, 5, 5, 5, 22, 22, 22, 21, 221, 211, 209, 214, 2143, 2097, 2099, 2101, 22236, 22342, 22185, 21901, 260380, 255719, 272745, 253801, 534993, 511377, 515131, 505316, 768498, 770864, 761531, 757728 ])
rnea_cpu_bf1 = rnea_cpu[0::4]

rnea_img = result_to_img(rnea_vtree, (8,4))

crba_cpu =np.array( [8, 8, 8, 7, 7, 256, 255, 141, 114, 79, 20137, 8846, 2369, 1672, 949, 526565, 107880, 30946, 23390, 16558, 2078497, 315392, 107734, 85784, 63488, 8301282, 863916, 324308, 272023, 222348, 32704701, 2407386, 1064404, 946915, 833427])
crba_cpu_bf1 = crba_cpu[0::5]
crba_cpu_bf1_01 = crba_cpu[1::5]
crba_cpu_bf1_1 = crba_cpu[2::5]
crba_cpu_bf1_2 = crba_cpu[3::5]
crba_cpu_bf2 = crba_cpu[4::5]

crba_vtree = np.array([2420, 2428, 2428, 2427, 2430, 2522, 2519, 2511, 2511, 2506, 2963, 2885, 2845, 2935, 2868, 7593, 5136, 4716, 4690, 4642, 21624, 8130, 6265, 6192, 6108, 78023, 15012, 10311, 10147, 9947, 14669508, 27922, 18919, 18434, 18056])
crba_gpu_bf1 = crba_vtree[0::5]
crba_gpu_bf1_01 = crba_vtree[1::5]
crba_gpu_bf1_1 = crba_vtree[2::5]
crba_gpu_bf1_2 = crba_vtree[3::5]
crba_gpu_bf2 = crba_vtree[4::5]


blocked_64 = np.array([1114, 1197, 2416, 3524, 21629, 41474, 62510])
work_efficient = np.array([ 625, 933, 2038, 8402, 64390, 131305, 260649])
normal_scan = np.array([937, 1163, 2393, 14639, 129984, 246428, 363029])

blocked_va_512 = np.array([ 206, 235, 362, 1430, 3862, 7224, 10808])

crba_img = result_to_img(crba_vtree, (7,5))

rnea_ns = np.flip(np.array([10, 100, 1000, 10000, 100000, 1000000, 2000000, 3000000]))
rnea_bfs = [1, 1.5 , 2, 1000]
crba_ns = np.flip(np.array([10, 100, 1000, 5000, 10000, 20000, 40000]))
crba_bfs = [1, 1.01, 1.1, 1.2, 2]
 
# print(crba_cpu_bf1, crba_gpu_bf1)
# barPlot(rnea_vtree_bf1, rnea_cpu_bf1, np.flip(rnea_ns), "Performance of RNEA")
barPlot(crba_gpu_bf1, crba_cpu_bf1, np.flip(crba_ns), "Performance of CRBA: Branching Factor = 1")
# barPlotForScans(blocked_64, work_efficient, normal_scan, np.flip(rnea_ns)[1:], "Performance of Scans")
# heatmap(rnea_img[:7], rnea_ns[:7], [1, 1.5 , 2, 1000], 'RNEA Heatmap', 0.96, 1.04)
# heatmap(crba_img, crba_ns, crba_bfs, 'CRBA Heatmap', 0.0, 1)
