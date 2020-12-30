import functools
import itertools
import numpy as np
from  math import pi

import time
import random

import sympy

import constants


def SYM_MOL_INIT(group_by,which_to_samp):
    # testing
    constants.DONE_ONCE = False
    constants.NUMBER_OF_DOWNSAMPLE = group_by
    constants.WHICH_DOWNSAMPLE = which_to_samp
    #constants.PRECOMPUTED_ANGLES = precomputed_angles

def CREATE_ANGLES_TEST(angoli, num_rot):
    # /////////////////////
    # COMBINATION WITH REPLACEMENT OF ANGLES

    pre_inspect = []
    for i in range(len(num_rot)):
        pre_inspect.append(angoli)
    #print("CEATE ANGLES")
    ret = list(itertools.product(*pre_inspect))
    ret = [list(x) for x in ret]
    #print(ret)
    return ret


def UPDATE_CONFORMATION(angles, conformation):
    #print("ANGLE NOW")
    #print(angles[conformation])
    return angles[conformation]




