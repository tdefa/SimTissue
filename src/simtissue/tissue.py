


from simtissue.arbitrary_shape import  generate_ellipse
import time
import random
from tqdm import tqdm
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from matplotlib import pyplot as plt
from scipy import sparse
from pathlib import Path
import sys
import tifffile
import numpy as np




def diffusion(distance_map, image_nuc, nucleus, speed=1):
    """ updade the distance map according to the initial object position and the speed
    """
    image = np.ones(image_nuc.shape)
    image[image_nuc != nucleus] = 0  # only zero expect at the nucleus area
    inverted_mask = np.ones(image.shape) - (image != 0).astype(int)
    distance_map = np.minimum(distance_map, ndi.distance_transform_edt(inverted_mask,
                                                                       sampling=speed))
    return distance_map




def simulate_single_cell_mask(mask_nuclei, cyto,
                              scale=np.array([3, 1.03, 1.03]),
                              proba_elipse=0,
                              intervals_speed=[[0.5, 0.8], [1.3, 4]],
                              rad_ellipse_range=[0.3, 1.5],
                              median_kernel=7,
                              random_seed = None):



    """
    apply individual distinction using parametrization
    :param masks_nuclei:
    :param cyto:
    :param scale:
    :param proba_elipse: probability of generating an elipse
    :param speed_range:
    :param rad_ellipse_range:
    :param median_kernel:
    :return:
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

        if random.uniform(0, 1) < proba_elipse:
            ellipse,  c, a , b = generate_ellipse(mask_nuclei, nuc,
                                                  rad_min=rad_ellipse_range[0],
                                                  rad_max=rad_ellipse_range[1])
            ellipse[mask_nuclei == nuc] = 0

            if len(np.unique(masks_nuclei[ellipse > 0])) != 1: ## if > 1 it intersect other nuclei
                ellipse = np.zeros(mask_nuclei.shape)
        else:
            ellipse = np.zeros(mask_nuclei.shape)
        distance_map = diffusion(distance_map, mask_nuclei + ellipse,
                                 nuc, speed=random_speed)


    labels = watershed(image=distance_map, markers=mask_nuclei, mask=cyto, compactness = 0)
    labels = ndi.median_filter(labels, size=median_kernel)
    return labels







