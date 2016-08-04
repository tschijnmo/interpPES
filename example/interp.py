#!/usr/bin/env python3

from interpPES import interp_PES
from ase.io import write

FORMAT = 'xyz'
INTERP = {
    'n_neb_images': 8,
    'n_interp_images': 0,
}
LINEAR_INTERP = {
    'n_neb_images': 8,
    'n_interp_images': 0,
    'interp': 'linear'
}
DUPL = 8

IMAGES = [
    {
        'file': 'reac1.xyz',
        'format': FORMAT, 'interp': LINEAR_INTERP, 'dupl': DUPL
    },
    {
        'file': 'reac3.xyz',
        'format': FORMAT, 'interp': INTERP, 'reorient': 1, 'dupl': DUPL
    },
    {
        'file': 'tss1.xyz',
        'format': FORMAT, 'interp': INTERP, 'reorient': 1, 'dupl': DUPL
    },
    {
        'file': 'prod.xyz',
        'format': FORMAT, 'reorient': 1, 'dupl': DUPL
    }
]

res = interp_PES(IMAGES)
write('res.pov', res, format='pov')
