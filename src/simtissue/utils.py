

import os
import time as time_module
from pathlib import Path

import numpy as np
import tifffile
from skimage.measure import regionprops
from tqdm import tqdm


def get_dict_coord_map(mask_cell,
                       ):
    dico_coord_map = {}
    list_cell = np.unique(mask_cell)
    if 0 in list_cell:
        list_cell = list_cell[1:]
    for cell in list_cell:
        out_ind = np.transpose(np.nonzero(mask_cell == cell))
        dico_coord_map[cell] = out_ind
    return dico_coord_map





def generate_dico_coord_map(cyto_path = "/home/tom/Bureau/phd/data3105_simulation/210426_repeat3/individual_cytoplasms_22_6_19_53/",
                            path_to_save = "/home/tom/Bureau/phd/data3105_simulation/210426_repeat3/dico_coord_map/",
                            extention = 'npy',
                            path_nuclei = '',
                            nuclei_mode = False):
    """
    generate dico of coordinate to then apply simulate_transcriptomic_expression
    :param cyto_path:
    :param path_to_save:
    :param path_nuclei:
    :param nuclei_mode: not implemented yet for nuclei
    :return:
    it saved dico_coord_map : { nuclei number : out_ind}
    """
    for folder_path in [path_to_save ]: # path_cyto_to_test is supposed to be already created
        try:
            os.mkdir(folder_path)
        except OSError as error:
            print(error)

    full_name_list = []
    for cyto in tqdm(list(Path(cyto_path).glob(f'*{extention}'))[:]):
        time_module.sleep(1)

        print(cyto)
        try:
            dico_coord_map = {}
            if 'npy' in extention:

                array_cyto = np.load(str(cyto))
            elif 'tif' in extention:
                array_cyto = tifffile.imread(str(cyto))
            else:
                raise(Exception(f'extention  {extention} not implemented'))
            list_nuclei = np.unique(array_cyto)
            if 0 in list_nuclei:
                list_nuclei = list_nuclei[1:]
            if len(list_nuclei) > 700:
                print('safe guard with 700')
                print(f"bug in the cyto generation ? {len(np.unique(cyto))}")
                continue
            if not nuclei_mode:
                for nuclei in list_nuclei:
                    print(nuclei)
                    out_ind = np.transpose(np.nonzero(array_cyto == nuclei))
                    dico_coord_map[nuclei] = out_ind
                np.save(path_to_save + str(cyto).split('/')[-1], dico_coord_map)
                full_name_list.append(str(cyto).split('/')[-1])
            else:
                raise(Exception('not implemented yet'))
        except Exception as e:
            print(e)
            print(f"error for file {cyto}")
    return full_name_list


def generate_dico_centroid(seg_mask):
    from skimage.measure import regionprops
    dico_centroid = {}
    props = regionprops(seg_mask)
    for pp in props:
        if len(pp.centroid) == 2:
            dico_centroid[pp.label] = [[0, pp.centroid[0], pp.centroid[1]]]
        else:
            assert len(pp.centroid) == 3
            dico_centroid[pp.label] = [pp.centroid]
    return dico_centroid



