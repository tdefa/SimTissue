

#%%
from tqdm import tqdm
import numpy as np
import random






############################################
# generate checker 2D board
############################################



############"
# file with fct to simulate various cytoplasmic simtissue
#############"


__all__ = ['checkerboard_mask', 'generate_ellipse', 'add_sphere_nuclei', "remove_nuclei"]

###### generate cheker-board mask

def checkerboard_mask(shape = [54, 1000, 1100],
                         cube_square_size = 100):
    """
    Generate checkerboard cell mask
    :param shape: dimension of the mask [z,y,x]
    :type shape: list[int]
    :param cube_square_size: size of the cube in pixel
    :return: checkerboard mask
    """


    assert len(shape) == 3
    m, n, d = shape[1], shape[2], cube_square_size
    checheckerboard = np.empty((m, n), dtype=int)
    arr_view = checheckerboard.reshape(m // d, d, n // d, d)
    vals = np.arange(m // d * n // d).reshape(m // d, 1, n // d, 1)
    arr_view[:] = vals
    checkerboard_mask = np.zeros(shape)
    for i in range(len(checkerboard_mask)):
        checkerboard_mask[i] = checheckerboard+1
    return  checkerboard_mask


def generate_ellipse(masks: np.ndarray,
                     cell_id : int,
                     radius_zyx :list  = None,
                     rad_min : float=0.5,
                     rad_max : float=2,
                     ):
    """
    generate nuclei as ellipse in the chosen cell
    :param masks: cell mask
    :type masks: np.ndarray
    :param cell_id: index of the chosen cell
    :type cell_id: int
    :param radius_zyx: list [Rz, Ry, Rx] of the fix parameter of the elipse leave it None to choose these parameter randomly between rad_min and rad_max
    :type raduis_zyx: list
    :param rad_min: Min value of the radius if radius_zyx is None
    :type rad_min: float
    :param rad_max:  Max value of the radius if radius_zyx is None
    :type rad_min: float

    :return: ellipse mask , ellipse parameter
    """

    if masks.ndim == 3 and len(masks) > 1:
        z = np.linspace(0, masks.shape[0], masks.shape[0])[:, None, None]
        y = np.linspace(0, masks.shape[1], masks.shape[1])[None, :, None]  # x values of interest
        x = np.linspace(0, masks.shape[2], masks.shape[2])[None, None, :]
        # y values of interest, as a "column" array
        list_cordo = list(zip(*np.nonzero((masks == cell_id).astype(int))))
        z0, y0, x0  = np.mean(list_cordo, axis=0)
        c,  b, a = (np.max(list_cordo, axis=0) - np.min(list_cordo, axis=0))

        if radius_zyx is None:
            c = c * random.uniform(rad_min, rad_max)
            a = a * random.uniform(rad_min, rad_max)
            b = b * random.uniform(rad_min, rad_max)
        else:
            c = c * radius_zyx[0]
            a = a * radius_zyx[1]
            b = b * radius_zyx[2]
        ellipse = ((z - z0) / c) ** 2 + ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 <= 1  # print((nuc, random_speed))
    elif masks.ndim == 3 and len(masks)==1:
        x = np.linspace(0, masks.shape[1], masks.shape[1])[None, :, None]  # x values of interest
        y = np.linspace(0, masks.shape[2], masks.shape[2])[None, None, :]
        # y values of interest, as a "column" array
        list_cordo = list(zip(*np.nonzero((masks == cell_id).astype(int))))
        z0, x0, y0 = np.mean(list_cordo, axis=0)
        c, a, b = (np.max(list_cordo, axis=0) - np.min(list_cordo, axis=0))
        if radius_zyx is None:
            a = a * random.uniform(rad_min, rad_max)
            b = b * random.uniform(rad_min, rad_max)
        else:
            a = a * radius_zyx[1]
            b = b * radius_zyx[2]
        ellipse = ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 <= 1  # print((nuc, random_speed))
    else:
        raise(Exception("Only implemented for 3D, for 2D use add a single z stack as [z,y,x]"))
    ellipse = ellipse * cell_id
    return ellipse, c, a, b


def add_sphere_nuclei(mask_cyto : np.ndarray,
                      cube_square_size :int=100,
                      nuclei_radius: int=25,
                      list_nuclei : list = []):

    """
    add nuclei as sphere in the chosen cell
    :param mask_cyto: cell mask
    :type mask_cyto: np.ndarray
    :param cube_square_size:   size of the cube in pixel
    :type cube_square_size: int
    :param nuclei_radius: radius of the nuclei in pixel in pourcentage of the cube_square_size
    :type nuclei_radius: float
    :param list_nuclei: list of nuclei to add
    :type list_nuclei: list
    :return: mask of the nuclei
    """

    rad = nuclei_radius / cube_square_size
    mask_nuclei = np.zeros(mask_cyto.shape)
    if list_nuclei is None:
        list_nuclei = np.unique(mask_cyto)
    for cell_id in list_nuclei:
        ellipse_nuc, c, a, b = generate_ellipse(
            masks=mask_cyto,
            cell_id=cell_id,
            rad_min=rad,
            rad_max=rad,
            radius_zyx=[rad, rad, rad]
        )
        mask_nuclei += ellipse_nuc
    return mask_nuclei


def remove_nuclei(
        mask_nuclei :np.ndarray,
        percent_to_remove : float = 0.5,
        list_nuc_to_keep : list = None):
    """
    remove nuclei from the existing nuclei mask to generate cell without nuclei like in smFISH experiment
    :param mask_nuclei: nuclei mask
    :type mask_nuclei: np.ndarray
    :param percent_to_remove: percent of nuclei to remove
    :type percent_to_remove: float
    :param list_nuc_to_keep: list of nuclei to keep
    :type list_nuc_to_keep: list
    :return: mask of the nuclei without some nuclei
    """
    list_nuc = np.unique(mask_nuclei)
    if 0 in list_nuc:
        assert list_nuc[0] == 0
        list_nuc = list_nuc[1:]
    from random import sample
    if list_nuc_to_keep is None:
        nb_nuc_to_keep = int(len(list_nuc) * (1 - percent_to_remove))
        list_nuc_to_keep = sample(list(list_nuc), nb_nuc_to_keep)
        assert len(list_nuc_to_keep) == nb_nuc_to_keep
    new_mask_nuclei = np.zeros(mask_nuclei.shape)
    for nuc in list_nuc_to_keep:
        new_mask_nuclei[mask_nuclei == nuc] = nuc
    return new_mask_nuclei







def elbow_grid_cube():
    """
    Function to generate grid of elbow shape. the size of the cell mask is hard coded to 1800x1600 and the size square cell is 100 pixel,
    the size of nuclei is 25 pixel
    :return: mask_cyto, mask_nuclei
    """

    ## simulate 2D and 3D grid
    shape = [1, 1800, 1600]
    step = 100
    mask_cyto = np.zeros(shape)

    m, n, d = shape[1], shape[2], step
    arr = np.empty((m, n), dtype=np.int64)
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
        #print(i)
        mask_nuclei[i] = arr + 1
    list_nuc = np.unique(mask_nuclei)
    #print(len(list_nuc))

    for nuc in tqdm(list_nuc[:]):
        #print(nuc)
        ellipse_nuc, c, a, b = generate_ellipse(
            masks = mask_nuclei,
            cell_id = nuc,
            rad_min=0.25,
            rad_max=0.25,
            radius_zyx=[0.25,0.25,0.25]

        )
        mask_nuclei_bis += ellipse_nuc
    #print(np.unique(mask_nuclei_bis))

    ## keep randomly one nuclei in every cytoplasm

    mask_nuclei_final = np.zeros(shape)
    list_unique_cyto = np.unique(mask_cyto)
    #print(len(list_unique_cyto))
    for cell in list_unique_cyto:
        candidate_nuc = np.unique(mask_nuclei_bis[mask_cyto==cell])[1:]
        #print(len(candidate_nuc))
        rint = random.randint(0, len(candidate_nuc)-1)

        mask_nuclei_final[mask_nuclei_bis == candidate_nuc[rint]] = cell


    return mask_cyto, mask_nuclei_final


