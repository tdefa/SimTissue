



import time
import os
import numpy as np
#from skimage.segmentation import watershed
#import tifffile
#from matplotlib import pyplot as plt
#from scipy import ndimage as ndi
import random
import argparse
from pathlib import Path
#from tqdm import tqdm





############################################
# generate checker 2D board
############################################



############"
# file with fct to simulate various cytoplasmic simtissue
#############"




###### generate cheker-board mask

def checkerboard_mask(shape = [54, 1000, 1100],
                         cube_square_size = 100):


    m, n, d = shape[1], shape[2], cube_square_size
    checheckerboard = np.empty((m, n), dtype=int)
    arr_view = checheckerboard.reshape(m // d, d, n // d, d)

    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    arr_view[:] = vals
    checkerboard_mask = np.zeros(shape)
    for i in range(len(checkerboard_mask)):
        checkerboard_mask[i] = checheckerboard+1
    return  checkerboard_mask


def generate_ellipse(masks,
                     nuc,
                     rad_min=0.5,
                     rad_max=2,
                     list_cab = None):
    """
    elipse for grid
    :param masks: 3D array
    :param nuc:
    :param rad_min:
    :param rad_max:
    :param list_cab:
    :return:
    """
    if masks.ndim == 3 and len(masks) > 1:
        z = np.linspace(0, masks.shape[0], masks.shape[0])[:, None, None]
        y = np.linspace(0, masks.shape[1], masks.shape[1])[None, :, None]  # x values of interest
        x = np.linspace(0, masks.shape[2], masks.shape[2])[None, None, :]
        # y values of interest, as a "column" array
        list_cordo = list(zip(*np.nonzero((masks == nuc).astype(int))))
        z0, y0, x0  = np.mean(list_cordo, axis=0)
        c,  b, a = (np.max(list_cordo, axis=0) - np.min(list_cordo, axis=0))
        if list_cab is None:
            c = c * random.uniform(rad_min, rad_max)
            a = a * random.uniform(rad_min, rad_max)
            b = b * random.uniform(rad_min, rad_max)
        else:
            c = c * list_cab[0]
            a = a * list_cab[1]
            b = b * list_cab[2]
        ellipse = ((z - z0) / c) ** 2 + ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 <= 1  # print((nuc, random_speed))
    elif masks.ndim == 3 and len(masks)==1:
        x = np.linspace(0, masks.shape[1], masks.shape[1])[None, :, None]  # x values of interest
        y = np.linspace(0, masks.shape[2], masks.shape[2])[None, None, :]
        # y values of interest, as a "column" array
        list_cordo = list(zip(*np.nonzero((masks == nuc).astype(int))))
        z0, x0, y0 = np.mean(list_cordo, axis=0)
        c, a, b = (np.max(list_cordo, axis=0) - np.min(list_cordo, axis=0))
        if list_cab is None:
            a = a * random.uniform(rad_min, rad_max)
            b = b * random.uniform(rad_min, rad_max)
        else:
            a = a * list_cab[1]
            b = b * list_cab[2]
        ellipse = ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 <= 1  # print((nuc, random_speed))
    else:
        raise(Exception("Only implemented for 3D, for 2D use add a single z stack as [z,y,x]"))
    ellipse = ellipse * nuc
    return ellipse, c, a, b


def add_sphere_nuclei(mask_cyto,
                      cube_square_size=100,
                      nuclei_radius=25,
                      list_nuclei=None):
    rad = nuclei_radius / cube_square_size

    mask_nuclei = np.zeros(mask_cyto.shape)
    if list_nuclei is None:
        list_nuclei = np.unique(mask_cyto)

    for nuc in list_nuclei:
        ellipse_nuc, c, a, b = generate_ellipse(
            masks=mask_cyto,
            nuc=nuc,
            rad_min=rad,
            rad_max=rad,
            list_cab=[rad, rad, rad]
        )
        mask_nuclei += ellipse_nuc
    return mask_nuclei


########## not use in tutorial

def generate_regular_mask(shape = [54, 1000, 1100],
                          cube_square_size = 100):

    """
    :param shape: [z,y,x]
    :param cube_square_size: cube_square_size
    :return:
    """
    # from tqdm import tqdm
    mask_nuclei = np.zeros(shape)
    mask_cyto = np.zeros(shape)
    mask_cyto_bis = np.zeros(shape)
    mask_nuclei_bis = np.zeros(shape)
    m, n, d = shape[1], shape[2], cube_square_size
    arr = np.empty((m, n), dtype=np.int)
    arr_view = arr.reshape(m // d, d, n // d, d)
    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    arr_view[:] = vals
    for i in range(len(mask_cyto)):
        mask_nuclei[i] = arr + 1
        mask_cyto[i] = arr + 1
    list_nuc = np.unique(mask_nuclei)
    print(len(list_nuc))
    for nuc in tqdm(list_nuc[:]):
        print(nuc)
        ellipse_nuc, c, a, b = generate_ellipse_for_regular_mask(masks = mask_nuclei,
                                                nuc = nuc,
                                                rad_min=0.25,
                                                rad_max=0.25,
                                                list_cab=[0.25,0.25, 0.25]
                                                                 )
        mask_nuclei_bis += ellipse_nuc
        ellipse_cyto, c, a, b = generate_ellipse_for_regular_mask(
                                                 mask_cyto,
                                                 nuc,
                                                 rad_min=0.5,
                                                 rad_max=0.5,
                                                 list_cab=[0.5,
                                                           0.5,
                                                           0.5]
                                                )
        mask_cyto_bis += ellipse_cyto
    return mask_nuclei, mask_cyto



def generate_ellipse_for_regular_mask(masks,
                                     nuc,
                                     rad_min=0.5,
                                     rad_max=2,
                                     list_cab = None):
    """
    elipse for grid
    :param masks: 3D array
    :param nuc:
    :param rad_min:
    :param rad_max:
    :param list_cab:
    :return:
    """
    if masks.ndim == 3 and len(masks) > 1:
        z = np.linspace(0, masks.shape[0], masks.shape[0])[:, None, None]
        y = np.linspace(0, masks.shape[1], masks.shape[1])[None, :, None]  # x values of interest
        x = np.linspace(0, masks.shape[2], masks.shape[2])[None, None, :]
        # y values of interest, as a "column" array
        list_cordo = list(zip(*np.nonzero((masks == nuc).astype(int))))
        z0, y0, x0  = np.mean(list_cordo, axis=0)
        c,  b, a = (np.max(list_cordo, axis=0) - np.min(list_cordo, axis=0))
        if list_cab is None:
            c = c * random.uniform(rad_min, rad_max)
            a = a * random.uniform(rad_min, rad_max)
            b = b * random.uniform(rad_min, rad_max)
        else:
            c = c * list_cab[0]
            a = a * list_cab[1]
            b = b * list_cab[2]
        ellipse = ((z - z0) / c) ** 2 + ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 <= 1  # print((nuc, random_speed))
    elif masks.ndim == 3 and len(masks)==1:
        x = np.linspace(0, masks.shape[1], masks.shape[1])[None, :, None]  # x values of interest
        y = np.linspace(0, masks.shape[2], masks.shape[2])[None, None, :]
        # y values of interest, as a "column" array
        list_cordo = list(zip(*np.nonzero((masks == nuc).astype(int))))
        z0, x0, y0 = np.mean(list_cordo, axis=0)
        c, a, b = (np.max(list_cordo, axis=0) - np.min(list_cordo, axis=0))
        if list_cab is None:
            a = a * random.uniform(rad_min, rad_max)
            b = b * random.uniform(rad_min, rad_max)
        else:
            a = a * list_cab[1]
            b = b * list_cab[2]
        ellipse = ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 <= 1  # print((nuc, random_speed))
    else:
        raise(Exception("Only implemented for 3D, for 2D use add a single z stack as [z,y,x]"))
    ellipse = ellipse * nuc
    return ellipse, c, a , b





def generate_regular_mask(shape = [54, 1000, 1100],
                          step= 100):

    """
    :param shape:
    :param step:
    :return:
    """
    from tqdm import tqdm
    mask_nuclei = np.zeros(shape)
    mask_cyto = np.zeros(shape)
    mask_cyto_bis = np.zeros(shape)
    mask_nuclei_bis = np.zeros(shape)
    m, n, d = shape[1], shape[2], step
    arr = np.empty((m, n), dtype=np.int)
    arr_view = arr.reshape(m // d, d, n // d, d)
    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    arr_view[:] = vals
    for i in range(len(mask_cyto)):
        mask_nuclei[i] = arr + 1
        mask_cyto[i] = arr + 1
    list_nuc = np.unique(mask_nuclei)
    print(len(list_nuc))
    for nuc in tqdm(list_nuc[:]):
        print(nuc)
        ellipse_nuc, c, a, b = generate_ellipse_for_regular_mask(masks = mask_nuclei,
                                                nuc = nuc,
                                                rad_min=0.25,
                                                rad_max=0.25,
                                                list_cab=[0.25,0.25, 0.25]
                                                                 )
        mask_nuclei_bis += ellipse_nuc
        ellipse_cyto, c, a, b = generate_ellipse_for_regular_mask(
                                                 mask_cyto,
                                                 nuc,
                                                 rad_min=0.5,
                                                 rad_max=0.5,
                                                 list_cab=[0.5,
                                                           0.5,
                                                           0.5]
                                                )
        mask_cyto_bis += ellipse_cyto
    return mask_nuclei, mask_cyto



#%%



palette = {"-1": "black",
           "0": "blue",
           "1": "orange",
           "2": "green",
           "3": "red",
           "4": "purple", "5": "brown", "6": "pink",
           "7": "olive",
           "8": "cyan",
           "9": "darkseagreen", "10": "coral", "11": 'lime', "12": "gold", "13": "lime",
           "14": "darkseagreen",
           "15": "chocolate", "16": "coral", "17": "plum", "18": "magenta",
           "19": "rosybrown", "20": "khaki", "21": "teal"}

for i in range(21, 500):
    palette[str(i)] = "#" + "%06x" % random.randint(0, 0xFFFFFF)

palette_int = {}
for kk in palette:
    palette_int[int(kk)] = palette[kk]




#%%

def elbow_grid():

    ## simulate 2D and 3D grid
    shape = [1, 1700, 1600]
    step = 100
    mask_cyto = np.zeros(shape)

    m, n, d = shape[1], shape[2], step
    arr = np.empty((m, n), dtype=np.int)
    #arr_view = arr.reshape(m // d, d, n // d, d)
    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    #arr_view[:] = vals
    mask_cyto[0][np.isin(arr, [0, 16, 17, 18])] = 1
    mask_cyto[0][np.isin(arr, [1, 2, 3, 3+16])] = 2
    mask_cyto[0][np.isin(arr, [4, 5, 6, 4+16])] = 3
    mask_cyto[0][np.isin(arr, [7, 5+16, 6+16, 7+16])] = 4
    mask_cyto[0][np.isin(arr, [32, 33, 34, 48])] = 5
    mask_cyto[0][np.isin(arr, [49, 50, 51, 35])] = 6
    mask_cyto[0][np.isin(arr, [36, 52, 53, 54])] = 7
    mask_cyto[0][np.isin(arr, [37, 38, 39, 55])] = 8
    for i in range(4):
        mask_cyto[0, 400 * i:400 * (i+1), 0:800 ] = mask_cyto[0, 0:400, 0:800 ] + 8*i
        mask_cyto[0, 400 * i:400 * (i+1), 800:1600 ] = mask_cyto[0, 0:400, 0:800 ] + 8*i + 32

    mask_nuclei = np.zeros(shape)
    mask_nuclei_bis = np.zeros(shape)

    for i in range(len(mask_nuclei)):
        print(i)
        mask_nuclei[i] = arr + 1
    list_nuc = np.unique(mask_nuclei)
    print(len(list_nuc))

    for nuc in tqdm(list_nuc[:]):
        print(nuc)
        ellipse_nuc, c, a, b = generate_ellipse_for_regular_mask(
            masks = mask_nuclei,
            nuc = nuc,
            rad_min=0.25,
            rad_max=0.25,
            list_cab=[0.25,0.25,0.25]
        )
        mask_nuclei_bis += ellipse_nuc

    ## keep randomly one nuclei in every cytoplasm

    mask_nuclei_final = np.zeros(shape)
    list_unique_cyto = np.unique(mask_cyto)
    print(len(list_unique_cyto))
    for cell in list_unique_cyto:
        candidate_nuc = np.unique(mask_nuclei_bis[mask_cyto==cell])[1:]
        print(len(candidate_nuc))
        rint = random.randint(0, len(candidate_nuc)-1)

        mask_nuclei_final[mask_nuclei_bis ==candidate_nuc[rint]] = cell

    import napari
    viewer = napari.Viewer()
    viewer.add_image(mask_cyto, name='mask')
    viewer.add_image(mask_nuclei_final, name='mask')

    viewer.add_labels(mask_cyto.astype(int), color=palette_int)  # , scale=[2, 8, 8]) # palette int from probe


    return mask_cyto, mask_nuclei_final




def elbow_grid_cube():

    ## simulate 2D and 3D grid
    shape = [1, 1800, 1600]
    step = 100
    mask_cyto = np.zeros(shape)

    m, n, d = shape[1], shape[2], step
    arr = np.empty((m, n), dtype=np.int)
    arr_view = arr.reshape(m // d, d, n // d, d)
    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    arr_view[:] = vals

    mask_cyto[0][np.isin(arr, [0, 16, 17, 18])] = 1
    mask_cyto[0][np.isin(arr, [1, 2, 3, 3+16])] = 2
    mask_cyto[0][np.isin(arr, [4, 5, 6, 4+16])] = 3
    mask_cyto[0][np.isin(arr, [7, 5+16, 6+16, 7+16])] = 4

    mask_cyto[0][np.isin(arr, [32])] = 5
    mask_cyto[0][np.isin(arr, [33])] = 6
    mask_cyto[0][np.isin(arr, [34])] = 7
    mask_cyto[0][np.isin(arr, [35])] = 8

    mask_cyto[0][np.isin(arr, [36])] = 9
    mask_cyto[0][np.isin(arr, [37])] = 10
    mask_cyto[0][np.isin(arr, [38])] = 11
    mask_cyto[0][np.isin(arr, [39])] = 12



    mask_cyto[0][np.isin(arr, [32+16, 33+16, 34+16, 48+16])] = 13
    mask_cyto[0][np.isin(arr, [49+16, 50+16, 51+16, 35+16])] = 14
    mask_cyto[0][np.isin(arr, [36+16, 52+16, 53+16, 54+16])] = 15
    mask_cyto[0][np.isin(arr, [37+16, 38+16, 39+16, 55+16])] = 16

    mask_cyto[0][np.isin(arr, [80])] = 17
    mask_cyto[0][np.isin(arr, [81])] = 18
    mask_cyto[0][np.isin(arr, [82])] = 19
    mask_cyto[0][np.isin(arr, [83])] = 20

    mask_cyto[0][np.isin(arr, [84])] = 21
    mask_cyto[0][np.isin(arr, [85])] = 22
    mask_cyto[0][np.isin(arr, [86])] = 23
    mask_cyto[0][np.isin(arr, [87])] = 24


    for i in range(3):
        mask_cyto[0, 600 * i:600 * (i+1), 0:800 ] = mask_cyto[0, 0:600, 0:800 ] + 24*i
        mask_cyto[0, 600 * i:600 * (i+1), 800:1600 ] = mask_cyto[0, 0:600, 0:800 ] + 24*i + 24 * 3

    mask_nuclei = np.zeros(shape)
    mask_nuclei_bis = np.zeros(shape)

    for i in range(len(mask_nuclei)):
        print(i)
        mask_nuclei[i] = arr + 1
    list_nuc = np.unique(mask_nuclei)
    print(len(list_nuc))

    for nuc in tqdm(list_nuc[:]):
        print(nuc)
        ellipse_nuc, c, a, b = generate_ellipse_for_regular_mask(
            masks = mask_nuclei,
            nuc = nuc,
            rad_min=0.25,
            rad_max=0.25,
            list_cab=[0.25,0.25,0.25]
        )
        mask_nuclei_bis += ellipse_nuc

    ## keep randomly one nuclei in every cytoplasm

    mask_nuclei_final = np.zeros(shape)
    list_unique_cyto = np.unique(mask_cyto)
    print(len(list_unique_cyto))
    for cell in list_unique_cyto:
        candidate_nuc = np.unique(mask_nuclei_bis[mask_cyto==cell])[1:]
        print(len(candidate_nuc))
        rint = random.randint(0, len(candidate_nuc)-1)

        mask_nuclei_final[mask_nuclei_bis ==candidate_nuc[rint]] = cell

    import napari
    if False:
            viewer = napari.Viewer()
            viewer.add_image(mask_cyto, name='mask')
            viewer.add_image(mask_nuclei_final, name='mask')

            viewer.add_labels(mask_cyto.astype(int), color=palette_int)  # , scale=[2, 8, 8]) # palette int from probe


    return mask_cyto, mask_nuclei_final




def elbow_grid_rectangle():

    ## simulate 2D and 3D grid
    shape = [1, 1800, 1600]
    step = 100
    mask_cyto = np.zeros(shape)

    m, n, d = shape[1], shape[2], step
    arr = np.empty((m, n), dtype=np.int)
    arr_view = arr.reshape(m // d, d, n // d, d)
    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    arr_view[:] = vals

    mask_cyto[0][np.isin(arr, [0, 16, 17, 18])] = 1
    mask_cyto[0][np.isin(arr, [1, 2, 3, 3+16])] = 2
    mask_cyto[0][np.isin(arr, [4, 5, 6, 4+16])] = 3
    mask_cyto[0][np.isin(arr, [7, 5+16, 6+16, 7+16])] = 4

    mask_cyto[0][np.isin(arr, [32, 33])] = 5
    mask_cyto[0][np.isin(arr, [34,35])] = 6

    mask_cyto[0][np.isin(arr, [36, 37])] = 7
    mask_cyto[0][np.isin(arr, [38, 39])] = 8



    mask_cyto[0][np.isin(arr, [32+16, 33+16, 34+16, 48+16])] = 9
    mask_cyto[0][np.isin(arr, [49+16, 50+16, 51+16, 35+16])] = 10
    mask_cyto[0][np.isin(arr, [36+16, 52+16, 53+16, 54+16])] = 11
    mask_cyto[0][np.isin(arr, [37+16, 38+16, 39+16, 55+16])] = 12

    mask_cyto[0][np.isin(arr, [80, 81])] = 13
    mask_cyto[0][np.isin(arr, [82, 83])] = 14

    mask_cyto[0][np.isin(arr, [84, 85])] = 15
    mask_cyto[0][np.isin(arr, [86, 87])] = 16


    for i in range(3):
        mask_cyto[0, 600 * i:600 * (i+1), 0:800 ] = mask_cyto[0, 0:600, 0:800 ] + 16*i
        mask_cyto[0, 600 * i:600 * (i+1), 800:1600 ] = mask_cyto[0, 0:600, 0:800 ] + 16*i + 16 * 3

    mask_nuclei = np.zeros(shape)
    mask_nuclei_bis = np.zeros(shape)

    for i in range(len(mask_nuclei)):
        print(i)
        mask_nuclei[i] = arr + 1
    list_nuc = np.unique(mask_nuclei)
    print(len(list_nuc))

    for nuc in tqdm(list_nuc[:]):
        print(nuc)
        ellipse_nuc, c, a, b = generate_ellipse_for_regular_mask(
            masks = mask_nuclei,
            nuc = nuc,
            rad_min=0.25,
            rad_max=0.25,
            list_cab=[0.25,0.25,0.25]
        )
        mask_nuclei_bis += ellipse_nuc

    ## keep randomly one nuclei in every cytoplasm

    mask_nuclei_final = np.zeros(shape)
    list_unique_cyto = np.unique(mask_cyto)
    print(len(list_unique_cyto))
    for cell in list_unique_cyto:
        candidate_nuc = np.unique(mask_nuclei_bis[mask_cyto==cell])[1:]
        print(len(candidate_nuc))
        rint = random.randint(0, len(candidate_nuc)-1)

        mask_nuclei_final[mask_nuclei_bis ==candidate_nuc[rint]] = cell

    import napari
    if False:
            viewer = napari.Viewer()
            viewer.add_image(mask_cyto, name='mask')
            viewer.add_image(mask_nuclei_final, name='mask')

            viewer.add_labels(mask_cyto.astype(int), color=palette_int)  # , scale=[2, 8, 8]) # palette int from probe


    return mask_cyto, mask_nuclei_final

















#### generate weak prior

def generate_weak_prior(mask_cyto,     min_range = -35,  max_range = 35):
    from scipy import ndimage as ndi
    from tqdm import tqdm
    import random


    unique_cyto = np.unique(mask_cyto)
    if 0 in unique_cyto:
        unique_cyto = unique_cyto[1:]
    random.shuffle(unique_cyto)
    mask_prior = np.zeros(mask_cyto.shape)

    for cell in tqdm(unique_cyto):
        kernel_size = random.randint(min_range, max_range)
        if kernel_size < 0:
            mask_prior[ndi.minimum_filter(mask_cyto == cell, size=-kernel_size)] = cell
        else:
            mask_prior[ndi.maximum_filter(mask_cyto == cell, size=kernel_size)] = cell

    return mask_prior

def folder_generate_week_prior(path_cyto = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/individual_cytoplasms_rd/",
                               path_save = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/individual_cytoplasms_rd_prior35_35/",
                               min_range = -35,
                               max_range = 35,
                               regex = "*tiff.npy",):

    Path(path_save).mkdir(parents=True, exist_ok=True)

    list_cyto = list(Path(path_cyto).glob(regex))
    list_cyto.reverse()
    for cyto in tqdm(list_cyto):
        print(cyto)
        mask_cyto = np.load(cyto)
        mask_prior = generate_weak_prior(mask_cyto, min_range, max_range)
        np.save(path_save+cyto.name[:-4], mask_prior)

    return None


#%%
if __name__ == "__main__":

    ## simulate 2D and 3D grid
    shape = [1, 1700, 1600]
    step = 100


    mask_cyto = np.zeros(shape)
    mask_nuclei = np.zeros(shape)
    mask_cyto_bis = np.zeros(shape)
    mask_nuclei_bis = np.zeros(shape)



    m, n, d = shape[1], shape[2], step
    arr = np.empty((m, n), dtype=np.int)
    arr_view = arr.reshape(m // d, d, n // d, d)
    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    arr_view[:] = vals




    for i in range(len(mask_cyto)):
        print(i)
        mask_nuclei[i] = arr + 1
        mask_cyto[i] = arr + 1

    list_nuc = np.unique(mask_nuclei)
    print(len(list_nuc))

    for nuc in tqdm(list_nuc[:]):
        print(nuc)
        ellipse_nuc, c, a, b = generate_ellipse_for_regular_mask(
            masks = mask_nuclei,
            nuc = nuc,
            rad_min=0.25,
            rad_max=0.25,
            list_cab=[0.25,0.25,0.25]
        )
        mask_nuclei_bis += ellipse_nuc

    #np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/cube2D_step100/cytoplasms/mask_cyto', mask_cyto)
    #np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/cube2D_step100/nuclei/mask_nuclei', mask_nuclei_bis)
    #np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/cube2D_step100/nuclei/mask_cyto', mask_nuclei_bis)


    #### elbow generation

    mask_cyto, mask_nuclei_final = elbow_grid()

    np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/elbow2D/ind_cytoplasms/mask_cyto0', mask_cyto)
    np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/elbow2D/nuclei/mask_nuclei0', mask_nuclei_final)

    for i in range(0,10):
        mask_cyto, mask_nuclei_final = elbow_grid()

        np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/elbow2D/ind_cytoplasms/mask_cyto' + str(i),
                mask_cyto)
        np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/elbow2D/nuclei/mask_cyto' + str(i),
                mask_nuclei_final)


    #### elbow generation 3Delbow_grid_rectangle

    for i in range(0, 10):
        mask_cyto, mask_nuclei_final = elbow_grid_cube()

        np.save('/media/tom/T7/regular_grid/simu1912/elbow_cube/ind_cytoplasms/mask_cyto' + str(i),
                mask_cyto)
        np.save('/media/tom/T7/regular_grid/simu1912/elbow_cube/nuclei/mask_cyto' + str(i),
                mask_nuclei_final)
    #### 3Delbow_grid_rectangle generation

    for i in range(0, 10):
        mask_cyto, mask_nuclei_final = elbow_grid_rectangle()

        np.save('/media/tom/T7/regular_grid/simu1912/elbow_rectangle/ind_cytoplasms/mask_cyto' + str(i),
                mask_cyto)
        np.save('/media/tom/T7/regular_grid/simu1912/elbow_rectangle/nuclei/mask_cyto' + str(i),
                mask_nuclei_final)


    import napari
    viewer = napari.Viewer()
    viewer.add_image(mask_cyto, name='mask')
    viewer.add_image(mask_nuclei_final, name='mask')

    viewer.add_labels(mask_cyto.astype(int), color=palette_int)



    #### weak prior generation
    from scipy import ndimage as ndi
    from tqdm import tqdm
    mask_cyto = np.load('/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/individual_cytoplasms_rd/01_CtrlNI_Chil3-Cy3_Mki67-Cy5_001.tiff.npy')

    min_range = -35
    max_range = 35



    unique_cyto = np.unique(mask_cyto)
    if 0 in unique_cyto:
        unique_cyto = unique_cyto[1:]
    random.shuffle(unique_cyto)
    mask_prior = np.zeros(mask_cyto.shape)

    for cell in tqdm(unique_cyto):
        kernel_size = random.randint(-min_range, max_range)
        if kernel_size < 0:
            mask_prior[ndi.minimum_filter(mask_cyto == cell, size=-kernel_size)] = cell
        else:
            mask_prior[ndi.maximum_filter(mask_cyto == cell, size=kernel_size)] = cell

    import napari

    viewer = napari.Viewer()

    viewer.add_image(mask_cyto, name='mask')
    viewer.add_image(mask_prior, name='mask_prior')


    #### elbow cube  generation 3D


    for i in range(0,10):
        mask_cyto, mask_nuclei_final = elbow_grid_cube()

        np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/elbow2D/ind_cytoplasms/mask_cyto' + str(i),
                mask_cyto)
        np.save('/home/tom/Bureau/phd/simulation/regular_grid/simu1912/elbow2D/nuclei/mask_cyto' + str(i),
                mask_nuclei_final)






