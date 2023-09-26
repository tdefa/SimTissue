




import anndata as ad
from scipy.sparse import csr_matrix
from .utils import get_dict_coord_map
import anndata as ad

from scipy import sparse



import os.path
import random
from scipy import ndimage
import argparse
import tifffile
from sklearn.utils.random import sample_without_replacement
import datetime
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import ssam
import itertools
from sklearn.neighbors import NearestNeighbors
from skimage import measure


import scanpy
from sklearn.utils.random import sample_without_replacement
from sklearn.utils import resample
import anndata as ad

from scipy import ndimage as ndi
from simtissue.utils import get_dict_coord_map



def create_fix_profile(dico_profile={'typeA': [120, 0],
                                     'typeB': [0, 10]},
                       cell_type_distribution={'typeA': list(range(0, 110, 2)),
                                               'typeB': list(range(1, 110, 2))}):
    dico_cell_type_label_rna = {}  # {cell_type : {cell number : expression vector}
    dico_cell_index = {}  # {cell number : {type:, rnaseq : expression vector}}
    for cell_type in cell_type_distribution:
        dico_cell_type_label_rna[cell_type] = {cell_id: dico_profile[cell_type] for cell_id in
                                               cell_type_distribution[cell_type]}
        for cell_id in cell_type_distribution[cell_type]:
            dico_cell_index[cell_id] = {'type': cell_type, "rnaseq": dico_profile[cell_type]}
    return dico_cell_type_label_rna, dico_cell_index

### simulate spots coordinate


def simulate_arbritrary_expression(dico_profile,
                  cell_type_distribution,
        genes_list_to_simulate = ['A', 'B'],
        image_name = "image0"
        ):

    dico_cell_type_label_rna, dico_cell_index = create_fix_profile(
        dico_profile=dico_profile,
        cell_type_distribution=cell_type_distribution
                )

    dict_coord_map = get_dict_coord_map(mask_cell=mask_cyto,
                                        )
    list_index_cell = list(dict_coord_map.keys())

    for cell in list_index_cell:
        dico_cell_index[cell]["ground_truth"] = {}
    list_list_gene = []
    list_list_coord = []
    for cell in list_index_cell:
        list_gene = []
        list_coord = []
        for gene in genes_list_to_simulate:
            gene_index = genes_list_to_simulate.index(gene)
            nb_rna = int(dico_cell_index[cell]["rnaseq"][gene_index])
            out_ind = dict_coord_map[cell]
            spots_postion = out_ind[sample_without_replacement(len(out_ind), min(nb_rna, len(out_ind)))]
            #dico_cell_index[cell]["ground_truth"][gene] = list(spots_postion)
            #dico_fish_channel[gene]["ground_truth"] += list(spots_postion)
            list_coord += list(spots_postion)
            list_gene += [gene] * len(spots_postion)
        list_list_gene.append(list_gene)
        list_list_coord.append(list_coord)
    ## gene expression vector

    list_expression_vector = []
    for cell_index in list_index_cell:
        list_expression_vector.append(dico_cell_index[cell_index]["rnaseq"])

    anndata = ad.AnnData(csr_matrix(list_expression_vector))
    anndata.var["features"] = genes_list_to_simulate
    anndata.var_names = genes_list_to_simulate
    anndata.obs["image_name"] = [image_name] * len(list_list_gene)
    anndata.obs["genes"] = list_list_gene
    anndata.obs["coordinate"] = list_list_coord
    anndata.obs["cell_index"] = list_index_cell
    return anndata


def filter_simulation(spots_position,
                      max_dist=3,
                      dict_scale={"x": 1, 'y': 1, "z": 1}):
    spots_position_scale = np.array(spots_position) * np.array([dict_scale['z'], dict_scale['y'], dict_scale["x"]])
    # print(len(spots_position_scale))
    nbrs = NearestNeighbors(n_neighbors=np.min([100, len(spots_position)]),
                            algorithm='ball_tree').fit(spots_position_scale)
    #ad = nbrs.kneighbors_graph(spots_position_scale)  ## can be optimize here
    distance = nbrs.kneighbors_graph(spots_position_scale, mode='distance')
    list_spots_too_take = []
    set_index_to_remove = set()
    for index_spot in range(len(spots_position_scale)):
        if index_spot not in set_index_to_remove:
            list_spots_too_take.append(index_spot)
            distance_index_spots = distance[index_spot].toarray()

            index_to_remove = np.nonzero(np.logical_and(distance_index_spots < max_dist, distance_index_spots != 0))[1]
            set_index_to_remove.update(set(list(index_to_remove)))

    return spots_position[list_spots_too_take]


def sim_spots_from_ref_anndata(

        ref_anndata,
        ind_cyto,
        selected_gene,
        image_name='',
        annotation_column="cell_ID",
        remove_neighbors=True,
        max_dist=0.3,
        dict_scale={"x": 0.103, 'y': 0.103, "z": 0.300},

):
    dico_coord_map = get_dict_coord_map(ind_cyto)
    random_indice = sample_without_replacement(len(ref_anndata),
                                               len(dico_coord_map))  # return a list
    list_index_cell = list(dico_coord_map.keys())

    scaling_factor = 3
    count_matrix = ref_anndata[random_indice, selected_gene].X
    print(ref_anndata[random_indice, selected_gene])
    list_cell_type = list(ref_anndata[random_indice, selected_gene].obs[annotation_column])

    if not isinstance(count_matrix, np.ndarray):
        count_matrix = count_matrix.toarray()

    list_list_gene = []
    list_list_coord = []
    list_expression_vector = []
    for cell_index in range(len(dico_coord_map)):
        cell = list(dico_coord_map.keys())[cell_index]
        # print(cell)
        ## get the expression vector
        expression_vector = count_matrix[cell_index] * scaling_factor
        simulated_expression_vector = []
        final_expression_vector = np.zeros(expression_vector.shape)
        list_gene = []
        list_coord = []
        for gene in selected_gene:
            nb_rna = expression_vector[list(selected_gene).index(gene)]

            out_ind = dico_coord_map[cell]
            spots_position = out_ind[sample_without_replacement(len(out_ind),
                                                                min(nb_rna, len(out_ind)))]
            if len(spots_position) > 0 and remove_neighbors:
                # print(len(spots_position))
                spots_position = filter_simulation(max_dist=3,
                                                   spots_position=spots_position,
                                                   dict_scale={"x": 1, 'y': 1, "z": 1})

            final_expression_vector[list(selected_gene).index(gene)] = len(spots_position)
            list_gene += [gene] * len(spots_position)
            list_coord += list(spots_position)

        list_list_gene.append(list_gene)
        list_list_coord.append(list_coord)
        list_expression_vector.append(final_expression_vector)

    anndata = ad.AnnData(sparse.csr_matrix(list_expression_vector))
    anndata.var["features"] = selected_gene
    anndata.var_names = selected_gene
    anndata.obs["image_name"] = [image_name] * len(list_list_gene)
    # anndata.obs["genes"] = np.array(list_list_gene)
    # anndata.obs["coordinate"] = np.array(list_list_coord)
    anndata.obs["cell_index"] = list_index_cell
    anndata.obs["cell_type"] = list_cell_type

    ### create csv file
    csv_list_z = []
    csv_list_y = []
    csv_list_x = []
    csv_list_gene = []
    csv_list_cell = []
    csv_list_cell_type = []
    for cell_index in range(len(list_index_cell)):
        if len(list_list_gene[cell_index]) > 0:
            csv_list_z += list(np.array(list_list_coord[cell_index])[:, 0])
            csv_list_y += list(np.array(list_list_coord[cell_index])[:, 1])
            csv_list_x += list(np.array(list_list_coord[cell_index])[:, 2])

            csv_list_gene += list_list_gene[cell_index]

            csv_list_cell += [list_index_cell[cell_index]] * len(list_list_gene[cell_index])

            csv_list_cell_type += [list_cell_type[cell_index]] * len(list_list_gene[cell_index])

    df_spots = pd.DataFrame({"z": csv_list_z,
                             "y": csv_list_y,
                             "x": csv_list_x,
                             "gene": csv_list_gene,
                             "cell": csv_list_cell,
                             'cell_type': csv_list_cell_type})

    anndata.uns["df_spots"] = df_spots

    return anndata, df_spots


if __name__ == '__main__':

    ### DEFINE THE TRANSCRIPTOMIC PROFILE OF EACH CELL
    cell_index_typeA = [i for i in range(1, 111) if i % 2 == 0]
    cell_index_typeB = [i for i in range(1, 111) if i % 2 == 0]
    genes_list_to_simulate = ['A', 'B']
    dico_cell_type_label_rna, dico_cell_index = create_fix_profile(
        dico_profile={'typeA': [50, 0],
                      'typeB': [0, 200]},
        cell_type_distribution={'typeA': cell_index_typeA,
                                'typeB': cell_index_typeB, }
    )

    ### simulate spots coordinate
    dict_coord_map = get_dict_coord_map(mask_cell=mask_cyto,
                                        )
    list_index_cell = list(dict_coord_map.keys())

    for cell in list_index_cell:
        dico_cell_index[cell]["ground_truth"] = {}
    list_list_gene = []
    list_list_coord = []
    for cell in list_index_cell:
        list_gene = []
        list_coord = []
        for gene in genes_list_to_simulate:
            gene_index = genes_list_to_simulate.index(gene)
            nb_rna = int(dico_cell_index[cell]["rnaseq"][gene_index] * dico_gene_scaling[gene])
            out_ind = dict_coord_map[cell]
            spots_postion = out_ind[sample_without_replacement(len(out_ind), min(nb_rna, len(out_ind)))]
            #dico_cell_index[cell]["ground_truth"][gene] = list(spots_postion)
            #dico_fish_channel[gene]["ground_truth"] += list(spots_postion)
            list_coord += list(spots_postion)
            list_gene += [gene] * len(spots_postion)
        list_list_gene.append(list_gene)
        list_list_coord.append(list_coord)

## gene expression vector
    list_expression_vector = []
    for cell_index in list_index_cell:
        list_expression_vector.append(dico_cell_index[cell_index]["rnaseq"])

    ### generate the list_expression_vectoranndata

    adata = ad.AnnData(csr_matrix(list_expression_vector))
    adata.var["features"] = genes_list_to_simulate
    adata.var_names = genes_list_to_simulate
    adata.obs["image_name"] = list_image_name
    adata.obs["genes"] = list_list_gene
    adata.obs["coordinate"] = list_list_coord
    adata.obs["cell_index"] = list_index_cell

