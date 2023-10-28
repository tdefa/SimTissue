

from pathlib import Path

import numpy as np
import scanpy
import pandas as pd

import tifffile
from anndata import AnnData
from tqdm import tqdm

#### select an  anndta IR5M and NI



path_save_ind_cyto = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/individual_cytoplasms_rd/"
path_to_save_anndata = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/NON_CONV_PAPER26/anndata/"
path_to_save_df = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/NON_CONV_PAPER26/dataframe/"
Path.mkdir(Path(path_to_save_anndata), exist_ok=True)
Path.mkdir(Path(path_to_save_df), exist_ok=True)

selected_gene = ['Atp6v0d2', 'Abcg1',# AM
             'Rtkn2',  'Igfbp2', #AT1
             'Sftpc','Cxcl15', #AT2,
            'Cd79a', #B_cells
             'Ms4a2', 'Fcer1a', #Basophils
             'Ccdc153', #Ciliated
             'Scgb3a2', 'Scgb1a1',#Club
             'Cst3',#DC
             'Cdh5', 'Clec14a',  #EC
             'Inmt', 'Pcolce2', # Fibroblasts
             'C1qc', 'C1qa', 'C1qb', # 'C3ar1', #IM
             'Upk3b',# Mesotheliocytes
             'Ifitm6','Plac8',# Monocytes
            'Ms4a4b', 'Ccl5', 'Hcst', # NK_T_cells
             'Gzma', 'Ncr1',# NK_cells
             'S100a9',# Neutrophils
             'Mmrn1',#Platelets
           'Acta2','Myh11', # SMC
             'Cd3g', 'Cd3d' #T_cells
             ]



from simtissue.control_expression import sim_spots_from_ref_anndata

# anndata = scanpy.read_h5ad( "/home/tom/Bureau/phd/markers_selection/data/20220321_lung_merge_27samples_raw_selected_with_subtype.h5ad")  # TODO
# anndata[anndata.obs['condition'] == 'IR_17Gy_5M', selected_gene].write_h5ad("/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/NON_CONV_PAPER26/IR5M.h5ad")
# anndata[np.isin(anndata.obs['condition'], [ 'NI_1M', 'NI_23M', 'NI_3M', 'NI_4M', 'NI_5M']), selected_gene].write_h5ad("/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/NON_CONV_PAPER26/NI.h5ad")

ref_anndata_NI = scanpy.read_h5ad("/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/NON_CONV_PAPER26/NI.h5ad")
ref_anndata_IR5M = scanpy.read_h5ad("/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/NON_CONV_PAPER26/IR5M.h5ad")


for path_ic in tqdm(list(Path(path_save_ind_cyto).glob("*.npy"))):
    print(path_ic.name)
    ind_cyto = np.load(path_ic)

    if 'NI' in path_ic.name:
        ref_anndata = ref_anndata_NI
    else:
        ref_anndata = ref_anndata_IR5M
    anndata, df_spots = sim_spots_from_ref_anndata(
        ref_anndata=ref_anndata,
        ind_cyto=ind_cyto,
        image_name=path_ic.stem,
        annotation_column="cell_ID",
        selected_gene=selected_gene,
        remove_neighbors=True,
        max_dist=0.3,
        dict_scale={"x": 0.103, 'y': 0.103, "z": 0.300},
    )
    anndata.write_h5ad(str(Path(path_to_save_anndata) / (path_ic.stem + '.h5ad')))
    df_spots.to_csv(Path(path_to_save_df) / (path_ic.stem + '.npy'))