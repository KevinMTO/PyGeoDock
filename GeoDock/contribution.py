
import sys
import constants


from mol2_nosymproc import *

from Rotation_Geo import rounding

import functools
import itertools
import numpy as np
from  math import pi

import time
import random




def get_contributes(M, distance_simplification):
    # CREATING SET OF ATOMS FIXED FROM WHICH WE WILL CALCULATE THE DISTANCE FROM
    # (atoms fixed with respect to center (in our case 5))
    molecules_sets = combo_influece_set(M)
    nodi = set(M.nodes())

    for ml in molecules_sets:
        subset = set(ml)
        nodi = nodi.difference(subset)

    SumDists = []
    ########################################################

    atoms_to_measure = molecules_sets
    atoms_to_measure.append(list(nodi))

    positions = M.nodes(data=True)  # RETURNS NODES WITH EXTRA DATA LIKE THE LABELS OF POSITIONS

    #print(positions)


    distances_list = {}

    # calcolo tutti i path con lunghezza maggiore di due
    # REMAKE SU IDs
    path_atoms_dict = {}
    iterate_via_ID = list(M.nodes)

    for IDx in iterate_via_ID:
        paths_from_IDx = paths_from(M, IDx)

        good_atom_for_measures = []
        for i in range(len(paths_from_IDx)):
            if (len(paths_from_IDx[i]) > 3):
                good_atom_for_measures.append(paths_from_IDx[i][len(paths_from_IDx[i]) - 1])
        good_atom_for_measures.sort()
        path_atoms_dict[IDx] = good_atom_for_measures




    # CALCULATING DISTANCES AND SUMMING UP THE CONTRIBUTIONS
    for a_set in range(len(atoms_to_measure) - 1):

        restante = a_set + 1

        for atomT in atoms_to_measure[a_set]:
            distances_list[atomT] = []

        for a_dist in range(restante, len(atoms_to_measure)):

            transformanda = atoms_to_measure[a_set]

            distance_from = atoms_to_measure[a_dist]

            #print("TRASFORMANDA")
            #print(transformanda)
            #print("distance_from")
            #print(distance_from)
################################ TWO ALTERNATIVES ##################################
            #######         ONE        ##############
            if distance_simplification == False:
                for atomT in transformanda:

                    totalcontribution = 0

                    for atomD in distance_from:
                        if (atomD in path_atoms_dict[atomT]):

                            #print("atomT prima di result")
                            #print(atomT)
                            result = distance_mapper(positions, atomT, atomD)
                            #print("result")
                            #print(result)
                            totalcontribution = totalcontribution + result
                            #print("totalcontribution")
                            #print(totalcontribution)
                            distances_list[atomT].append(atomD)

                    SumDists.append((totalcontribution))

            #######         TWO        ##############
            """
            if distance_simplification == True:

                # check che ditance_from e trasformanda non siano vuoti
                if (len(transformanda) > 0 and len(distance_from) > 0):

                    # for atomT in random.sample(transformanda, min(3,len(transformanda))):

                    # totalcontribution = 0

                    # for atomD in random.sample(distance_from, min(3,len(transformanda))):

                    # considero solo l'atomo che sta a metà lista (sia per atomT che atomD)
                    atomT = transformanda[int(len(transformanda) / 2)]

                    totalcontribution = 0

                    atomD = distance_from[int(len(distance_from) / 2)]

                    if (atomD in path_atoms_dict[atomT]):
                        result = distance_mapper(positions, atomT, atomD)

                        totalcontribution = totalcontribution + result

                        distances_list[atomT].append(atomD)

                    SumDists.append((totalcontribution))
            """

################################ SUMMING CONTRIBUTIONS #########################
    #print("SumDists")
    #print(SumDists)
    #print(type(SumDists))

    contributi_sommati = functools.reduce(lambda x, y: x + y, SumDists)


    #print("contributi sommati")
    #print(contributi_sommati)

    return contributi_sommati


