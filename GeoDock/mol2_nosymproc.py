
from biopandas.mol2 import PandasMol2

import pandas as pd


from io import StringIO
import csv

from Rotation_Geo import *

from networkx import *

import functools

import random
import constants

from operator import itemgetter
#############################################################################################################################
#IS TESING PHASE UP?



SAMPLED = []

def testing_dowsnsampling(cut_l):
    if(constants.NUMBER_OF_DOWNSAMPLE<=constants.SIZE_OF_SAMPLE and constants.TEST_PHASE_MD==True):

        clnew = list(itemgetter(*constants.WHICH_DOWNSAMPLE)(cut_l))

        #CASO IN CUI ITEMGETTER PRENDE UN SOLO ELEMENTO
        if(type(clnew[0])!=tuple):
            clnew=[(clnew[0],clnew[1])]

        if(constants.DONE_ONCE==False):
            constants.DONE_ONCE = True
            global SAMPLED
            SAMPLED = clnew
            constants.MAP_ang_rot = list(zip( SAMPLED, constants.ANGLES_NOW ))
            #print("constants.MAP_ang_roT")
            #print(constants.MAP_ang_rot)


        return SAMPLED
    elif(constants.TEST_PHASE_MD==False):
        return cut_l
    else:
        raise StopIteration
#############################################################################################################################


def mol2_atoms(file):
    parsed_mol = PandasMol2().read_mol2(file).df
    parsed_mol.drop(parsed_mol[parsed_mol['atom_name'].str.contains('H')].index, inplace=True)
    return parsed_mol

def require_bond(ID ,bonds):
    """ returns tuple of atoms of bond """
    b = bonds[bonds['bond_id' ]==ID]
    b = b[['origin_atom_id', 'target_atom_id']]
    return [tuple(x) for x in b.to_numpy()]



def get_coord(data ,ID):
    """returns coordinates of atom by ID"""
    return ((data[data['atom_id' ]==ID][['x' ,'y' ,'z']].to_numpy()).tolist())[0]  # ho tolto /100!!!!



def mol2_bonds(filename, atoms_frame):
    """Returns dataframe of bonds of mol2 file"""
    with open(filename, 'r') as f_input:
        begin = False
        rows = []
        export = []
        for line in f_input:
            if(line.startswith("@<TRIPOS>BOND") and begin==False):
                begin = True
            elif(line.startswith("@<TRIPOS>") and begin==True):
                break
            elif begin:
                rows.append(next(csv.reader(StringIO(line))))
        for row in rows:
            export.append(row[0].split())

    ##print(export)
    for dfline in export:
        for i in range((5-len(dfline))):
            dfline.append('NaN')
    ##print(export)
    df = pd.DataFrame(export, columns = ['bond_id', 'origin_atom_id', 'target_atom_id', 'bond_type','STATUS'])

    df = df.astype({'bond_id' :'int64' ,'origin_atom_id' :'int64' ,'target_atom_id': 'int64', 'bond_type': 'U','STATUS':'U'})

    df.drop(df.loc[~df['target_atom_id'].isin(atoms_frame['atom_id'])].index, inplace=True)
    df.drop(df.loc[~df['origin_atom_id'].isin(atoms_frame['atom_id'])].index, inplace=True)

    return df



#############################################################################################################################

def mol_as_graph(pmol, bonds):
    """Returns Graph from networkx with verteces and edges coherent with mol2 file. Verteces contain their position as
        label
    """
    M = Graph()

    for row in pmol.itertuples():
        M.add_node(row.atom_id)


    edges = []
    edge_type = {}

    for row in bonds.itertuples():
        couple = (row.origin_atom_id, row.target_atom_id)

        edges.append(couple)
        edge_type[couple] = row.bond_type

    M.add_edges_from(edges)
    set_edge_attributes(M, edge_type, 'bond_type')

    nodes = list(M.nodes)
    positions = dict()

    for a in nodes:
        poly_pos = atom_PolyM(pmol, a)
        positions[a] = poly_pos

    set_node_attributes(M, positions, 'coords')

    return M


############################################################################################################

def rotables_cleansing(M, rotables):
    result = []
    for i in range(len(rotables)):
        if(M.edges[rotables[i][0] ,rotables[i][1]]['bond_type'] in ['1', 'ar']):

            result.append(rotables[i])

    return result

def get_rotables(M):
    """Returns rotable bonds in graph. Found by removing and checking if connected iteratively."""
    if (constants.DONE_ONCE == True):
        return SAMPLED

    cuts = []
    for edge in M.edges():
        G = Graph(M)
        G.remove_edge(edge[0], edge[1])
        if (not (is_connected(G))):
            set1 = node_connected_component(G, edge[0])
            set2 = node_connected_component(G, edge[1])
            if (len(set1) != 1 and len(set2) != 1):
                cuts.append(edge)
    cuts.sort()
    cuts = rotables_cleansing(M,cuts)
    #=========================================================
    #ADDED FOR TESTING PURPOSES BUT WE CAN REMOVE IT

    constants.SIZE_OF_SAMPLE  = len(cuts)
    #######################################
    
    #cuts_new li ordina in base alla centrality
    atoms_btw=atoms_btw_ctrl(M)
    atoms_btw_ctrl_list=[atoms_btw[i][0] for i in range(len(atoms_btw))]
    
    cuts_new=[]
    for i in atoms_btw_ctrl_list:
        for j in cuts:
            if ( i in j):
                if (j in cuts_new):
                    pass
                else:
                    cuts_new.append(j)

    ####cuts = testing_dowsnsampling(cuts)
    cuts = testing_dowsnsampling(cuts_new)

    #=========================================================
    return cuts


def atoms_btw_ctrl(G):
    """Function to get atoms betweeness centrality"""
    sort_nodes_bc = sorted(betweenness_centrality(G).items(), key=lambda x: x[1], reverse=True)
    return sort_nodes_bc


def paths_from(G, n):
    """Returns shortest path from the atoms with greatest betweeness centrality"""
    PATHS_from_o = [v for v in shortest_path(G, source=n).values()]
    return PATHS_from_o


def routes_by_rotable(M):
    """Returns which which combination of rotables is implicit in the shortest path to an atom."""
    n = atoms_btw_ctrl(M)[0][0]

    list_PATHS_from_o = paths_from(M, n)

    cuts = get_rotables(M)

    linked_rots = []
    for i in list_PATHS_from_o:

        combs = []

        for j in cuts:
            if (j[0] in i and j[1] in i):
                if (i.index(j[0]) > i.index(j[1])):
                    combs.append((j[1], j[0]))
                else:
                    combs.append((j[0], j[1]))
        linked_rots.append(combs)

    return linked_rots


def all_cuts_combo(M):
    """Returns all legal possible combinations of rotables"""

    linked_rots = routes_by_rotable(M)

    linked_rots = [el for el in linked_rots if (len(el))]  # alcuni shortest path non contengono rotabili

    for x in linked_rots:
        x.reverse()

    cleanlist = []
    for x in linked_rots:
        if (x not in cleanlist):
            cleanlist.append(x)

    return cleanlist  #####################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OCCHIOOOOOOOOOOO
    # return [[(5, 4)], [(4, 2), (5, 4)]] #CASO SEMPLIFICATO FORZATO


#############################################################################################################################

def all_rots(M,pmol):
    """Returns a list of the matrices identifying the legal possible combinations of rotables.
       Sorted in the same way as the output of 'all_cuts_combo'"""

    cleanlist = all_cuts_combo(M)

    matrici = []

    solo_bonds=[x[0] for x in constants.MAP_ang_rot]


    for combined in cleanlist:
        R = []
        for legame in combined:
            
            ##SYMBOLIC CALCULATIONS ARE LOST

            #TODO FOR NOW IS AN EASY FIX WITH GLOBAL VARIABLES BLAH

            #print("legame for matrice")
            #print(legame)
            esito = [(legame[0] in y and legame[1] in y)for y in solo_bonds]
            index = esito.index(True)
            angle_inserted = constants.MAP_ang_rot[index][1]
            R.append(matrix( angle_inserted,
                             get_coord(pmol, legame[0]), get_coord(pmol, legame[1])))


            #print("MATRICE")
            #print(R)

        matrici.append(functools.reduce(lambda r1, r2: r1.dot(r2), R))

    return matrici


def combo_influece_set(M):
    """Returns list of sets of atoms divided by which cuts combination influeces their postion.
       Sorted in the same order of all_cuts_combo"""
    transformanda = []

    C = atoms_btw_ctrl(M)[0][0]

    paths = paths_from(M, C)

    routes = routes_by_rotable(M)

    for x in routes:
        x.reverse()

    combos = all_cuts_combo(M)

    for index in range(len(combos)):
        transformanda.append([])

    for j in range(len(routes)):
        for index in range(len(combos)):
            if (routes[j] == combos[index]):
                transformanda[index].append(paths[j][-1])

    # #print(combos)
    for index in range(len(combos)):
        # #print(len(combos[index]))
        if (not (not (combos[index][0][1] in transformanda[index]))):

            if (len(combos[index]) > 1):
                transformanda[index].append(combos[index][0][1])

            elif (combos[index][0][1] in transformanda[index]):
                transformanda[index].remove(combos[index][0][1])

    results = []
    for x in transformanda:
        atomset = []
        for atom in x:
            if (not (atom in atomset)):
                atomset.append(atom)

        results.append(atomset)

    return results  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # return [[2], [1, 3,14,15,16]] #CASO SEMPLIFICATO FORZATO


#############################################################################################################################

def atom_PolyM(pmol, atom):
    """Returns Vector of Polynomials representing the position of an atom"""
    atomcoords = get_coord(pmol, atom)
    atomcoords.append(1)
    atomcoords = numpy.asarray(atomcoords)
    return atomcoords


def move_atoms(M,pmol):
    """Modifies the position of all atoms through the application of the right Matrix on the right set of atoms"""
    atoms_to_move = combo_influece_set(M)
    rots_cuts = all_cuts_combo(M)
    R_list = all_rots(M, pmol)

    positions = M.nodes(data=True)

    for composed_track in range(0, len(rots_cuts)):

        transformanda = atoms_to_move[composed_track]

        for atomT in transformanda:
            positions[atomT]['coords'] = R_list[composed_track].dot(positions[atomT]['coords'])



###########################################################
#### ATTENTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###########################################################


def get_SymDistance(I, atomT): #USED FOR PARALLELIZING
    """Returns difference vector of two atoms"""
    # builder = [ i*i for i in (I - atomT)]
    # Dist = PolyMatrix(builder)
    # return builder
    #print("I")
    #print(I)
    #print("atomT")
    #print(atomT)
    return (I - atomT)


def distance_mapper(positions, atomTcoords, I_coord):
    """Returns the distance of two atoms identified by their atom number"""
    Dist = get_SymDistance(positions[I_coord]['coords'], positions[atomTcoords]['coords'])
    #print("Dist")
    #print(Dist)
    result = ( (Dist[0]*Dist[0])+(Dist[1]*Dist[1])+(Dist[2]*Dist[2]) )**(1/2)

    return result
