








import anndata as ad
from scipy.sparse import csr_matrix
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




from tqdm import tqdm



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
    for index_spot in tqdm(range(len(spots_position_scale))):
        if index_spot not in set_index_to_remove:
            list_spots_too_take.append(index_spot)
            distance_index_spots = distance[index_spot].toarray()

            index_to_remove = np.nonzero(np.logical_and(distance_index_spots < max_dist, distance_index_spots != 0))[1]
            set_index_to_remove.update(set(list(index_to_remove)))

    return spots_position[list_spots_too_take] , list_spots_too_take, set_index_to_remove








if __name__ == '__main__':



    dico_commu = np.load("/home/tom/Bureau/phd/simulation/mycode/dico_folder/dico_dico_commulist_cubeABpruning1cube2D_step100dico_simulation_A100_B100_AB200.npy",
            allow_pickle=True).item()
    path_df_folder = "/media/tom/T7/regular_grid/simu1912/cube2D_step100/dataframe_folder/A100_B100_AB200/"
    path_df_folder_save = "/media/tom/T7/regular_grid/simu1912/cube2D_step100/dataframe_folder/A100_B100_AB200_max3/"


    for image_name in tqdm(dico_commu):
        df_spots_label = dico_commu[image_name]['df_spots_label']

        df_spots_label.to_csv(path_df_folder + image_name[ : -4] +'.csv', index=False)


    ##### save df from dico




    for path_df in tqdm(list(Path(path_df_folder).glob('*'))):


        df_spots = pd.read_csv(path_df)

        selected_gene = np.unique(df_spots['gene'].values)

        total_list_x = []
        total_list_y = []
        total_list_z = []
        total_list_gene = []
        total_list_cell = []
        total_list_celltype = []

        for gene in selected_gene:
            list_x = np.array(df_spots[df_spots['gene'] == gene].x)
            list_y = np.array(df_spots[df_spots['gene'] == gene].y)
            list_z = np.array(df_spots[df_spots['gene'] == gene].z)
            list_gene = np.array(df_spots[df_spots['gene'] == gene].gene)
            list_cell = np.array(df_spots[df_spots['gene'] == gene].cell)
            list_celltype = np.array(df_spots[df_spots['gene'] == gene].cell_type)

            spots_position = df_spots[df_spots['gene'] == gene][['z', 'y', 'x']].values


            spots_position, list_spots_too_take, index_to_remove = filter_simulation(
                                                spots_position,
                                                  max_dist=3,
                                                  dict_scale={"x": 1, 'y': 1, "z": 1})


            total_list_x += list(list_x[list_spots_too_take])
            total_list_y += list(list_y[list_spots_too_take])
            total_list_z += list(list_z[list_spots_too_take])
            total_list_gene += list(list_gene[list_spots_too_take])
            total_list_cell += list(list_cell[list_spots_too_take])
            total_list_celltype += list(list_celltype[list_spots_too_take])



            print(f'len index_to_remove {len(index_to_remove)} for gene {gene}')

        df_spots = pd.DataFrame({'x': total_list_x,
                                    'y': total_list_y,
                                    'z': total_list_z,
                                    'gene': total_list_gene,
                                    'cell': total_list_cell,
                                    'cell_type': total_list_celltype})

        df_spots.to_csv(path_df_folder_save + path_df.name, index=False)

        print()

    ############## from df folder to dico_dico_comm ##############


    path_df_folder_save = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/dataframe_cluster_filter_max2/"

    dico_dico_commu = {}
    for path_df in tqdm(list(Path(path_df_folder_save).glob('*'))):
        df_spots = pd.read_csv(path_df)
        dico_dico_commu[path_df.stem + ".tiff.npy"] = {}
        dico_dico_commu[path_df.stem+ ".tiff.npy"]["df_spots_label"] = df_spots


    np.save(path_df_folder_save + "dico_dico_commu_filter_max2", dico_dico_commu)






