
import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage as ndi
import skimage


def plot_contour(mask_cyto=None,
                 mask_nuclei=None,
                 figsize=(15, 15),
                 anndata=None,
                 dico_color={"A": "r", "B": "b"},
                 spot_size=3,
                 linewidth=1):
    fig, ax = plt.subplots(figsize=figsize)

    if mask_cyto is not None:
        mask_cyto = np.amax(mask_cyto, axis=0)
        cyto_list = np.unique(mask_cyto)
        for nuc_id in cyto_list:
            countour = skimage.measure.find_contours(mask_cyto, nuc_id)
            try:
                plt.plot(countour[0][:, 1], countour[0][:, 0], linewidth=linewidth, c='k')
            except Exception as e:
                print(e)
                print("no contour for cell", nuc_id)

    if mask_nuclei is not None:
        mask_nuclei = np.amax(mask_nuclei, axis=0)
        cyto_list = np.unique(mask_nuclei)
        for nuc_id in cyto_list:
            countour = skimage.measure.find_contours(mask_nuclei, nuc_id)
            try:
                plt.plot(countour[0][:, 1], countour[0][:, 0], linewidth=linewidth, c='k')
            except Exception as e:
                print(e)
                print("no contour for cell", nuc_id)

    if anndata is not None:
        list_list_gene = list(anndata.obs["genes"])
        list_list_coord = list(anndata.obs["coordinate"])
        list_unique_gene = np.unique(np.concatenate(list_list_gene))
        dico_gene_coord = {}
        for gene in list_unique_gene:
            dico_gene_coord[gene] = []
        for gene, coord in zip(np.concatenate(list_list_gene),
                               np.concatenate(list_list_coord)):
            dico_gene_coord[gene].append(coord)

        for gene, arr in dico_gene_coord.items():
            arr = np.array(arr)
            ax.scatter(arr[:, 2], arr[:, 1],
                       c=dico_color[gene],
                       # cmap=list_cmap1,,  #  gene_color_dico[gene],   #'#%02X%02X%02X' % (r(),r(),r()),
                       s=spot_size)
    return fig, ax
