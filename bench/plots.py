import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib as mpl

# Taken from: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
def heatmap(xs, yaxis, xaxis, title, vmin, vmax):

    fig, ax = plt.subplots()
    im = ax.imshow(xs, vmin=vmin, vmax=vmax)

    ax.set_xlabel('Branching Factor')
    ax.set_ylabel('Number Of Bodies')

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, )
    cbar.set_ticks(ticks=[vmin, 1, vmax], labels=[vmin, 1, vmax])
    cbar.ax.set_ylabel(r'Row isolated deviation from BF = 1 ($\frac{bf_i}{bf_1}$)', rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(range(len(xaxis)), labels=xaxis,
                  rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(yaxis)), labels=yaxis)

    # Loop over data dimensions and create text annotations.
    for i in range(len(yaxis)):
        for j in range(len(xaxis)):
            text = ax.text(j, i, xs[i, j],
                           ha="center", va="center", color="w")

    ax.set_title(title)
    fig.tight_layout()
    plt.show()


rnea_res = [1263, 1261, 1266, 1273, 1973, 1963, 1963, 1955, 2900, 2922, 2917, 2919, 6602, 6630, 6614, 6629, 32313, 32427, 32374, 32358, 257080, 259655, 258470, 257950]
rnea_img = np.resize(np.array(rnea_res), (6,4))
rnea_img2 = []
for x in rnea_img:
    bf1_res = x[0]
    rnea_img2.append([round(a/bf1_res, 4) for a in x])
rnea_img2 = np.flip(np.array(rnea_img2), 0)
rnea_bf1 = rnea_res[0::4]
rnea_bf2 = rnea_res[1::4]
rnea_bf10 = rnea_res[2::4]
rnea_bf1000 = rnea_res[3::4]
crba_res = [1462, 1448, 1452, 1449, 2271, 2260, 2248, 2264, 3420, 3361, 3365, 3347, 9294, 6315, 6352, 6338, 26444, 8868, 8835, 8919, 92460, 14866, 14814, 14822, 246129, 29544, 29430, 29321]
crba_bf1 = crba_res[0::4]
crba_bf2 = crba_res[1::4]
crba_bf10 = crba_res[2::4]
crba_bf1000 = crba_res[3::4]

rnea_ns = np.flip(np.array([10, 100, 1000, 10000, 100000, 1000000]))
crba_ns = np.flip(np.array([10, 100, 1000, 5000, 10000, 20000, 40000]))
 
heatmap(rnea_img2, rnea_ns, [1, 2, 10, 1000], 'RNEA heatmap', 0.98, 1.02)
