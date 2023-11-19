

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
    :param shape: [z,y,x]
    :type shape: list[int]
    :param cube_square_size:
    :return:
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
    :param radius_zyx: list [Rz, Ry, Rx] of the fix parameter of the elipse leave it None to choose these parameter randomly
    :type raduis_zyx: list
    :param rad_min: if radius_zyx is None min value of the radius
    :type rad_min: float
    :param rad_max: if radius_zyx is None max value of the radius
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
    :param nuclei_radius:
    :param list_nuclei:
    :return:
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
    remove nuclei from the existing nuclei segmentation mask to generate cell without nuclei like in smFISH experiment
    :param mask_nuclei:
    :param percent_to_remove:
    :param list_nuc_to_keep:
    :return:
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






