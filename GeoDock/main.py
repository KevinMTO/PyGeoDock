import builtins
import sys

import faulthandler
faulthandler.enable()

import argparse
import constants


from mol2_nosymproc import *

from contribution import *
from inject_solution import *

import functools
import itertools
import numpy
from  math import pi

import time
import random


if __name__ == "__main__":

    text = 'lol non abbiamo una descrizione'

    # Create the parser and add arguments
    parser = argparse.ArgumentParser(description=text)

    parser.add_argument(dest='is_test', choices=['True', 'False'], help="activate 'a scaglioni' ")

    parser.add_argument(dest='filepath', type=str, help="mol2 input file path")
    parser.add_argument(dest='logname', type=str, help="logname")
    parser.add_argument(dest='iteration', type=int, help="iteration for experiment")

    parser.add_argument(dest='group_by', type=int, help="how many")


    parser.add_argument(dest='distance_simplification', choices=['True', 'False'], help="flag for loop perforation")
    parser.add_argument(dest='multiple_steps', choices=['True', 'False'], help="flag for 'a scaglioni' ")
    parser.add_argument(dest='toggle_brute_hubo', choices=['True', 'False'], help="flag for choosing hubo or inspection: True is by inspection")

    parser.add_argument(dest='gran', type=int, help="integer by which we divide 180 degrees")
    parser.add_argument(dest='round_value', type=int, help="rounding values")


    args = parser.parse_args()

    print(args)

    #================================================================================================================================
    #================================================================================================================================
    group_by = args.group_by
    iteration = args.iteration

    is_test = (args.is_test == 'True')

    filepath = args.filepath
    logname = args.logname

    distance_simplification = (args.distance_simplification == 'True')
    #se distance_simplification=True allora considero solo un atomo per  ogni corpo rigido


    # se step>0 e multiple_steps=True,  fissa i rotabili del run precedente e ottimizza solo i successivi
    MULTIPLE_STEPS = (args.multiple_steps == 'True')
    constants.MULTIPLE_STEPS=MULTIPLE_STEPS

    toggle_brute_hubo= (args.toggle_brute_hubo == 'True')

    d = args.gran

    #decidi approx degli angoli
    round_val = args.round_value


############################################## START CODE ############################################################


    # =============================================================================
    #                                   SET-UP
    # =============================================================================


    constants.TEST_PHASE_MD = False


    SYM_MOL_INIT(0, ())# create and initialize global variables with default values

    log = logname+'#s'+str(group_by)+'#it'+str(iteration)+'.txt'
    print(log)


    ##########################################


    ############################################
    #           SETTINGS FOR ACCURACY

    # devi settare quanti angoli vuoi (4=45, 6=30,18=10) dividi l'intervallo in pi/d parti
    #d = 2
    angoli = [pi * k / d for k in range(1, 2*d+1 )] #range(1, 2*d+1)


    print("angles= ")
    print(str(angoli))
    print(" \n")
    print(" \n")

    for i in range(len(angoli)):
        angoli[i]=round(angoli[i],round_val)

    ####################################
    # PARSE THE INITIAL MOLECULE

    pmol = mol2_atoms(filepath)

    my_bonds = mol2_bonds(filepath,pmol)

    print("molecule \n")
    print(pmol.to_string())
    print(" \n")
    print(my_bonds.to_string())
    print(" \n")
    print(" \n")


    ####################################
    # BUILD MOLECULE GRAPH
    M = mol_as_graph(pmol,my_bonds)

    ####################################
    ROTABLES = get_rotables(M)
    

    if(group_by > len(ROTABLES)):
        print("Invalid group by- greater than num rotables")
        sys.exit(-1)
        
    ####################################
    #  DIVIDING THE BONDS IN SMALL GROUPS TO PROCESS IN BATCHES
    ROTABLES_SECTIONING=[]
    counter_modulus=0
    group = []
    for num in range(len(ROTABLES)):
        group.append(num)
        counter_modulus += 1

        if(counter_modulus>=group_by):
            ROTABLES_SECTIONING.append(group)
            group=[]
            counter_modulus=0
    if(len(group)!=0):
        ROTABLES_SECTIONING.append(group)
        
    ####################################
    print("ROTABLES_SECTIONING")
    print(ROTABLES_SECTIONING)

    
    
    
    
    constants.TEST_PHASE_MD = is_test #TODO se possibile rimuovere ormai inutile
    
    

    # =============================================================================

    #                                   CONTRIBUTI

    # =============================================================================
    start_time_global = time.time()
    angoli_risultato=[]

    ##MAXIFYING VARIABLES
    max_val = 0

    for i in range(3):
        for group_of_bonds in ROTABLES_SECTIONING:

            constants.ALL_ANGLES = CREATE_ANGLES_TEST(angoli, group_of_bonds)
            #print("ALL_ANGLES")
            #print(constants.ALL_ANGLES)

            ##MAXIFYING VARIABLES
            max_rot = []
            for z in range(len(group_of_bonds)):
                max_rot.append(0.0)

            start_time_local = time.time()
            for conformation in range(len(constants.ALL_ANGLES)):
                T = nx.Graph(M)

                constants.ANGLES_NOW = UPDATE_CONFORMATION(constants.ALL_ANGLES, conformation)
                #LOGGING
                #print("ANGLES_NOW")
                #print(constants.ANGLES_NOW)


                try:


                    sampling = len(group_of_bonds)
                    #print("SAMPLING")
                    #print(sampling)

                    SYM_MOL_INIT(sampling, group_of_bonds)

                    move_atoms(T, pmol)

                    #print("posizioni variate")
                    #print(T.nodes(data=True))

                except StopIteration:
                    print("finished experiment iteration")
                else:

                    # =============================================================================
                        #GeoDock
                    # =============================================================================

                    val = get_contributes(T, distance_simplification )
                    #valuto le combinazioni e trovo quella che da il massimo e me la salvo

                    #print("val: "+ str(val)+ " vs "+str(max_val))# LOGGING

                    if(val > max_val):
                        max_rot = constants.ANGLES_NOW
                        max_val = val

                    #print("contributions %s seconds ---" % (time.time() - start_time_local))
                    #start_time_local = time.time()


            print("--- inspection %s seconds ---" % (time.time() - start_time_local))


            constants.ANGLES_NOW = max_rot
            move_atoms(M, pmol) #trasformo il grafo secondo la configurazione che mi da il massimo
            angoli_risultato.append(max_rot)
            print("ANGLES: "+str(max_rot) +" MAX_VAL: " +str(max_val))
        angoli_risultato.append([-1])



    print("--- GLOBALLY %s seconds ---" % (time.time() - start_time_global))
    print(" BEST ANGLES")
    print(angoli_risultato)
