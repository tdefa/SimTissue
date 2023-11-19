
import numpy as np
from matplotlib import pyplot as plt
import skimage


def plot_contour(mask_cyto=None,
                 mask_nuclei=None,
                 figsize=(15, 15),
                 anndata=None,
                 dico_color={"A": "r", "B": "b"},
                 spot_size=3,
                 linewidth=1):
    fig, ax = plt.subplots(figsize=figsize)

    if mask_cyto.ndim == 3:
        mask_cyto = np.amax(mask_cyto, axis=0)
    if mask_nuclei.ndim == 3:
        mask_nuclei = np.amax(mask_nuclei, axis=0)


    if mask_cyto is not None:
        cyto_list = np.unique(mask_cyto)
        for nuc_id in cyto_list:
            countour = skimage.measure.find_contours(mask_cyto, nuc_id)
            try:
                plt.plot(countour[0][:, 1], countour[0][:, 0], linewidth=linewidth, c='k')
            except Exception as e:
                print(e)
                print("no contour for cell", nuc_id)

    if mask_nuclei is not None:
        from scipy import ndimage as ndi
        contour_nuclei = (mask_nuclei  > 0).astype(int) - ndi.minimum_filter((mask_nuclei > 0).astype(int), size=3)
        contour_nuclei = np.array(list(zip(*np.nonzero(contour_nuclei))))
        plt.scatter(contour_nuclei[:, 1], contour_nuclei[:, 0],
                    s=linewidth,
                    linewidths = linewidth, c='black')

    ### plot cell border
    contour_nuclei = [[0, i] for i in
                      range(mask_nuclei.shape[-1])] + [[i, 0] for i
                                                  in range(mask_nuclei.shape[-2])] + [[mask_nuclei.shape[-2] - 1, i] for i
                                                                                 in range(mask_nuclei.shape[-1])] + [
                         [i, mask_nuclei.shape[-1] - 1] for i
                         in range(mask_nuclei.shape[-2])]
    contour_nuclei = np.array(contour_nuclei)
    plt.scatter(contour_nuclei[:, 1], contour_nuclei[:, 0],
                linewidths=linewidth,
                s=linewidth,
                c='black')

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
