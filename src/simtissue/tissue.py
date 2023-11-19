


import time
import random
from tqdm import tqdm
from scipy import ndimage as ndi
from skimage.segmentation import watershed
import numpy as np


__all__ = ['simulate_single_cell_mask']

def diffusion(distance_map, image_nuc, nucleus, speed=1):
    """ updade the distance map according to the initial object position and the speed
    """
    image = np.ones(image_nuc.shape)
    image[image_nuc != nucleus] = 0  # only zero expect at the nucleus area
    inverted_mask = np.ones(image.shape) - (image != 0).astype(int)
    distance_map = np.minimum(distance_map, ndi.distance_transform_edt(inverted_mask,
                                                                       sampling=speed))
    return distance_map




def simulate_single_cell_mask(mask_nuclei : np.ndarray,
                              cyto: np.ndarray = None,
                              scale :  np.ndarray = np.array([3, 1.03, 1.03]),
                              intervals_speed: list[list] = [[0.5, 0.8], [1.3, 4]],
                              median_kernel=7,
                              random_seed = None):



    """
    generate a single cell mask from a nuclei.  Individual cytoplasms are defined by growing cells from segmented nuclei.
     Each cell grows at random speed to add irregularity in the cell size.
    :param masks_nuclei: nuclei segmentation mask
    :type masks_nuclei: np.ndarray
    :param cyto: (Optional), cytoplasm segmentation mask to simulate area without cell
    :type cyto: np.ndarray
    :param scale :  scale of the image in z,y,x
    :type scale :  np.ndarray
    :param speed_range: speed range of the cell growth
    :type speed_range: list
    :param median_kernel: median kernel size to smooth the cell mask
    :return: cell mask
    """

    if random_seed is None:
        random_seed = int(str(time.time())[-7:].split('.')[-1])
        random.seed(random_seed)

    nuclei_list = np.unique(mask_nuclei)[1:]
    distance_map = np.ones(mask_nuclei.shape) * np.inf
    for nuc in tqdm(list(nuclei_list)[:]):
        t = time.time()

        random_speed = np.array([random.uniform
            (*random.choices(intervals_speed, weights=[r[1 ] -r[0] for r in intervals_speed])[0]) * scale[i] for i in range(len(scale))])
        distance_map = diffusion(distance_map, mask_nuclei,
                                 nuc, speed=random_speed)
    labels = watershed(image=distance_map, markers=mask_nuclei, mask=cyto, compactness = 0)
    cell_mask = ndi.median_filter(labels, size=median_kernel)
    return cell_mask







