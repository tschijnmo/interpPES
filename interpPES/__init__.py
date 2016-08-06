"""
Utilities for interpolating between stationary points on PES
============================================================

"""

import itertools

from ase.build import minimize_rotation_and_translation
from ase.neb import NEB
from ase.optimize import MDMin
from ase.io import read


#
# Interpolation by NEB method
# ---------------------------
#


def interp_by_neb(initial, final, n_images, interp='idpp',
                  calculator=None, constraint=None, fmax=0.5, steps=10,
                  **kwargs):
    """Interpolate between the initial and final points by NEB method.

    :param initial: The initial image.
    :param final: The final image.
    :param int n_images: The number of image between the initial and final
        images.
    :param string interp: The interpolation method.

    :param Callable calculator: The callable to generate calculators for the
        images.  It can be set to None if no real optimization is intended.
    :param constraint: The constraint for the images.
    :param fmax: The maximal force to stop the optimization.
    :param steps: The number of optimization steps allowed.

    :param kwargs:  All other keyword arguments are forwarded to the
        initializer of the ASE NEB class.

    :return: The list of images between the points.
    """

    # To circumvent the error from interpolation when we have no middle images.
    if n_images == 0:
        return [initial, final]

    images = [initial]
    for _ in range(n_images):
        image = initial.copy()
        images.append(image)
        if calculator is not None:
            image.set_calculator(calculator())
        if constraint is not None:
            image.set_constraint(constraint)
        continue
    images.append(final)

    neb = NEB(images, **kwargs)
    neb.interpolate(method=interp)

    if calculator is not None:
        dyn = MDMin(neb)
        dyn.run(fmax=fmax, steps=steps)

    return neb.images


def interp_by_multiscale_neb(initial, final, n_neb_images, n_interp_images,
                             **kwargs):
    """Interpolate by multi-scale NEB method.

    By multi-scale NEB method, it is meant that first a full NEB interpolation
    is performed between the images, then pure IDPP interpolation is performed
    between the computed NEB images.
    """

    neb_res = interp_by_neb(initial, final, n_neb_images, **kwargs)

    return itertools.chain.from_iterable(
        interp_by_neb(i, j, n_interp_images)
        for i, j in zip(neb_res, neb_res[1:])
    )


#
# Stationary point argumentation
# ------------------------------
#


def aug_point(point, n_dupl):
    """Augment a stationary point by making it a list of n duplications."""

    return (point for _ in range(n_dupl))


#
# The main driver
# ---------------
#


def interp_PES(images, flatten=True):
    """Interpolate between an series of images.

    The images should be specified by a dictionary, with keys,

    file

        The name of the file to read the image from.

    format

        The format of the file for the image, according to the requirement of
        the ``read`` function in the ASE IO module.

    atoms

        The ASE atoms object for the image point.  If it is present, the file
        will not be attempted to be read.

    dupl

        The number of duplication of the image.

    interp

        The dictionary of arguments to the function
        :py:func:`interp_by_multiscale_neb`.  For each of the images, it gives
        the method of interpolation between the current image and its
        successor.  No interpolation is performed in the absence of it.

    reorient

        If the image point should be reorientated with respect to the previous
        image point.

    """

    FILE = 'file'
    FORMAT = 'format'
    ATOMS = 'atoms'
    DUPL = 'dupl'
    INTERP = 'interp'
    REORIENT = 'reorient'

    # To hold the duplication of each image as well as the path to its
    # successor.
    FRAMES = '_frames'

    # Read the images if it is needed.
    imgs = []
    for i in images:
        img = dict(i)
        if ATOMS in i:
            atoms = i[ATOMS]
        else:
            if FORMAT not in i:
                atoms = read(i[FILE])
            else:
                atoms = read(i[FILE], format=i[FORMAT])
        if REORIENT in i:
            minimize_rotation_and_translation(imgs[-1][ATOMS], atoms)
        img[ATOMS] = atoms
        img[FRAMES] = [atoms, ]
        imgs.append(img)
        continue

    # Make duplication for the images.
    for i in imgs:
        if DUPL in i:
            i[FRAMES].extend(aug_point(i[ATOMS], i[DUPL]))
        continue

    # Interpolate between the stationary points.
    for curr, succ in zip(imgs, imgs[1:]):
        if INTERP not in curr:
            continue
        else:
            curr[FRAMES].extend(interp_by_multiscale_neb(
                curr[ATOMS], succ[ATOMS], **curr[INTERP]
            ))

    frame_sets = (i[FRAMES] for i in imgs)
    if flatten:
        return list(itertools.chain.from_iterable(frame_sets))
    else:
        return list(frame_sets)
