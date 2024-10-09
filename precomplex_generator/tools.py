#!/usr/bin/env python
#
# This file is part of the precomplex_generator
#
# Copyright (C) BASF SE 2022
#
# The precomplex_generator is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The precomplex_generator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the precomplex_generator. If not, see <https://www.gnu.org/licenses/>.

from __future__ import division, print_function, unicode_literals
try:
    from pymatgen.core import Molecule
    from pymatgen.io.babel import BabelMolAdaptor
    # from StructureGeneration.symmetry_tools import *
except:
    pass

from openbabel import openbabel as ob
import os, sys
import networkx as nx
import numpy as np
import itertools
import fnmatch
import traceback
from element_dict import element_dict

# from openbabel import OBElementTable -> OBElementTable now removed
# // OB 2.x
# OBElementTable etab;
# const char *elem = etab.GetSymbol(6);
# unsigned int atomic_num = etab.GetAtomicNum(elem);

# // OB 3.0
# #include <openbabel/elements.h>
# const char *elem = OBElements::GetSymbol(6);
# unsigned int atomic_num = OBElements::GetAtomicNum(elem);


def ShortestPathSearch(obmol, a1, a2):

    """
    Search for shortest path between two atoms (a1 and a2) in an openbabel molecule (obmol).
    """

    edges = []
    for bond in ob.OBMolBondIter(obmol):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetIdx() > atom2.GetIdx() and (atom2.GetIdx(), atom1.GetIdx()) not in edges:
            edges.append((atom2.GetIdx(), atom1.GetIdx()))
        elif atom2.GetIdx() > atom1.GetIdx() and (atom1.GetIdx(), atom2.GetIdx()) not in edges:
            edges.append((atom1.GetIdx(), atom2.GetIdx()))

    G = nx.Graph()
    G.add_edges_from(edges)

    path = nx.shortest_path(G, a1, a2)

    return len(path), path

def RingSearch(obmols, atIDs):

    """
    Search for rings in 1 or 2 openbabel molecules
    obmol[0] has atID1 and atID2, obmol[1] (if it exists) has atID3 and atID4
    """

    # intramolecular reaction
    if len(obmols) == 1:
        obmol1 = obmols[0]["obm"]
        obmol2 = obmols[0]["obm"]
        atID1 = atIDs[0][0]
        atID2 = atIDs[0][1]
        len_path1, path1 = ShortestPathSearch(obmol=obmol1, a1=atID1, a2=atID2)
        ring_size = len_path1
        rings = [path1]
    else:
        # bimolecular reaction
        obmol1 = obmols[0]["obm"]
        obmol2 = obmols[1]["obm"]
        atID1 = atIDs[0][0]
        atID2 = atIDs[0][1]
        atID3 = atIDs[1][0]
        atID4 = atIDs[1][1]

        #print("now evaluate ShortestPath")
        len_path1, path1 = ShortestPathSearch(obmol=obmol1, a1=atID1, a2=atID2)
        len_path2, path2 = ShortestPathSearch(obmol=obmol2, a1=atID3, a2=atID4)

        ring_size = len_path1 + len_path2
        rings = [path1, path2]

    return ring_size, rings


def center_of_points(array=None):

    """
    Takes a (nx3) numpy-array as input and returns the center of points.
    If only a single point is given it is returned.

    :param arr: Numpy array with the structure coordinates.
    :return: Center of points (cop).
    """

    length = len(array)

    if length > 0:
        x = ((np.sum(array[:, 0])) / length)
        y = ((np.sum(array[:, 1])) / length)
        z = ((np.sum(array[:, 2])) / length)
        cop = np.asarray([x, y, z])
    else:
        cop = array

    return cop


def angle_func(v1, v2):

    """
    Takes two vectors (v1 and v2) and returns the angle between them in degrees.
    """

    a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    a = np.asscalar(a)
    if a > 1:
        b = 0
    elif a < -1:
        b = np.pi
    else:
        b = np.degrees(np.arccos(a))

    return b


def norm_vec(v):

    """
    Takes a vector (v) as input and returns it normed to unit length.
    """

    a = np.linalg.norm(v)
    if a > 0:
        b = v / a
    else:
        b=v

    return b


def rot_mat(v1, v2):

    """
    Takes two vectors as input and returns the rotation matrix to rotate vector 1 onto vector 2.

    :param v1: Vector that is to be rotated.
    :param v2: Vector to which the alignment should lead.
    :return: r: (3x3) Rotation Matrix as numpy array.
    """

    v = np.cross(v1, v2)
    c = np.dot(v1, v2)

    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])
    if abs(1+c) > 1e-3:
        r = np.identity(3) + vx + (np.matmul(vx, vx) * (1 / (1 + c)))
    else:
        r = -np.identity(3)

    return r


def dihedral_func(x):

    """
    Takes a (4x3) matrix containing the coordinates of the four considered points as input and returns their dihedral
    angle in degrees.

    :param p: Numpy-array with the coordinates of the four points.
    :return: Dihedral angle in degrees.
    """

    x0, x1, x2, x3 = x[0:4]

    v0 = -1.0 * (x1 - x0)
    v1 = x2 - x1
    v2 = x3 - x2

    # normalize v1 so that it does not influence magnitude of vector
    # rejections that come next
    v1 /= np.linalg.norm(v1)

    # vector rejections
    # v = projection of v0 onto plane perpendicular to v1
    #   = v0 minus component that aligns with v1
    # w = projection of v2 onto plane perpendicular to v1
    #   = v2 minus component that aligns with v1
    v = v0 - np.dot(v0, v1) * v1
    w = v2 - np.dot(v2, v1) * v1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(v1, v), w)

    return np.degrees(np.arctan2(y, x))


def kov_rad(symbol1, symbol2):

    """
    Determine "covalent radius" between two atoms (plus 0.3 Angstrom)
    """

    rad = {"H": 0.32, "He": 0.32,
           "Li": 0.76, "Be": 0.45, "B": 0.82, "C": 0.77, "N": 0.71, "O": 0.73, "F": 0.71, "Ne": 0.69,
           "Na": 1.02, "Mg": 0.72, "Al": 0.535, "Si": 1.11, "P": 1.06, "S": 1.02, "Cl": 0.99, "Ar": 0.97,
           "K": 1.38, "Ca": 1.00, "Sc": 1.44, "Ti": 1.36, "V": 1.25, "Cr": 1.27, "Mn": 1.39, "Fe": 1.25,
           "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31, "Ga": 1.26, "Ge": 1.22, "As": 1.21, "Se": 1.16,
           "Br": 1.14, "Kr": 1.10,
           "Rb": 1.66, "Sr": 1.32, "Y": 1.62, "Zr": 1.48, "Nb": 1.37, "Mo": 1.45, "Tc": 1.31, "Ru": 1.26,
           "Rh": 1.35, "Pd": 1.31, "Ag": 1.53, "Cd": 1.48, "In": 1.44, "Sn": 1.41, "Sb": 1.38, "Te": 1.35,
           "I": 1.33, "Xe": 1.30,
           "Cs": 2.25, "Ba": 1.98, "La": 1.69, "Hf": 1.50, "Ta": 1.38, "W": 1.46, "Re": 1.59, "Os": 1.28,
           "Ir": 1.37, "Pt": 1.38, "Au": 1.36, "Hg": 1.49, "Tl": 1.48, "Pb": 1.46, "Bi": 1.46, "Po": 1.40,
           "At": 1.45, "Rn": 1.45}

    #Radii of elements found at http://de.wikipedia.org/wiki/Kovalenter_Radius
    #Ionic radii for alkaline and alkaline earth metals + Al found in ref: https://doi.org/10.1107/S0567739476001551

    rad_1 = rad[symbol1]
    rad_2 = rad[symbol2]

    #ToDo: check factor
    #rad_out = 1.3*(rad_1+rad_2)
    rad_out = rad_1+rad_2+0.3

    return rad_out

# really needed???
def kov_radius(asym):

    """
    Takes the atomic symbol of the element as input and returns the atomic radius. The value for hydrogen is doubled,
    for the alkaline and alkaline earth metals and Al, the ionic radii are used. The numbers are taken from:
    -> For Al increased by 50%, otherwise bonds not identified

    http://de.wikipedia.org/wiki/Kovalenter_Radius
    https://doi.org/10.1107/S0567739476001551
    """

    kovr_data = {'H': 0.64, 'Li': 0.76, 'Be': 0.45, 'B': 0.82, 'C': 0.77, 'N': 0.71, 'O': 0.73, 'F': 0.71, 'Na': 1.02,
                 'Mg': 0.72, 'Al': 0.803, 'Si': 1.11, 'P': 1.06, 'S': 1.02, 'Cl': 0.99, 'K': 1.38, 'Ca': 1.00,
                 'Sc': 1.44, 'Ti': 1.36, 'V': 1.25, 'Cr': 1.27, 'Mn': 1.39, 'Fe': 1.25, 'Co': 1.26, 'Ni': 1.21,
                 'Cu': 1.38, 'Zn': 1.31, 'Ga': 1.26, 'Ge': 1.22, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Rb': 1.66,
                 'Sr': 1.32, 'Y': 1.62, 'Zr': 1.48, 'Nb': 1.37, 'Mo': 1.45, 'Tc': 1.31, 'Ru': 1.26, 'Rh': 1.35,
                 'Pd': 1.31, 'Ag': 1.53, 'Cd': 1.48, 'In': 1.44, 'Sn': 1.41, 'Sb': 1.38, 'Te': 1.35, 'I': 1.33,
                 'La': 1.80, 'Bi': 1.46}

    # doubled for H atom

    radius = kovr_data[asym]

    return radius

def file_read_list(filename=None):

    """
    Takes a file, opens it and returns its contents as a list of lists of strings.
    If the input is a xyz-File, the first two lines are removed, if it is a tors.dat-File from a Fortran-Code,
    the keywords $tors and $end are removed.

    :param filename: File to be read into program, including extension.
    :return: List of list of strings of the file content.
    """

    with open(filename, 'r') as fileOpened:
        lines = fileOpened.readlines()
        list_out = []
        for line in lines:
            list_out.append(line.split())

    if filename[-8:] == 'tors.dat':
        list_out = list_out[1:-1]

    elif filename[-4:] == '.xyz':
        list_out = list_out[2:]

    return list_out


def file_read_dict(filename=None):

    """
    Takes a file, opens it and returns its contents as a dictionary of lists of strings.

    :param filename: File to be read into program, including extension.
    :return: Dictionary of list of strings of the file content.
    """

    try:
        with open(filename, 'r') as fileOpened:
            lines = fileOpened.readlines()
            dict_out = {}
            for line in lines:
                line_split = line.split()
                dict_out[line_split[0]] = list(map(int, [int(x) for x in line_split[3:]]))

    except:
        dict_out = {}

    return dict_out


def input_gen(new_data=None, data=None, case=None):

    """
    Takes values for the first/second add or break indices or bonds to be rotated as new_data and adds the respective
    keywords depending on the case for the input into the Fortran-Code. Rotation step width should be set in a more
    differentiated way so that

    :param new_data: New input (atom IDs as integers or strings).
    :param data: Data is input and will be returned with the new data appended.
    :param case: Case to determine the added keywords: fadd, add, breaks, tors.
    :return: Data is returned with the relevant ids and keywords added to it.

    Input_gen and file_write so seemingly complicated to cope with all possible cases.
    Should be replaced in the future by much simpler direct argument passing to Fortran subroutines...
    """

    if case == 'fadd':
        data.append(['$fadd'])

    elif case == 'add':
        data.append(['$add'])

    elif case == 'breaks':
        if ['$breaks'] not in data:
            data.append(['$breaks'])

    elif case == 'tors':
        if ['$tors'] not in data:
            data.append(['$tors'])

    if len(data) == 1:
        if len(new_data) == 3:
            data.append(new_data)
        elif len(new_data) > 3:
            for single_data in new_data:
                data.append(single_data)
        else:
            print('Error in data appending after first keyword: ', new_data)

    else:
        if case == 'fadd':
            data.append(new_data)

        elif case == 'add':
            data.append(new_data)

        elif case == 'breaks':
            for single_break in new_data:
                if len(single_break) > 0:
                    data.append(single_break)

        elif case == 'tors':

            # More sophisticated handling of rotatable bonds to remove redundancies

            old_temp = []
            to_dels = []

            for old in data:
                if len(old) == 3:
                    old_temp.append(sorted(old[0:2]))

            new_temp = []
            for new in new_data:
                new_temp.append(sorted(list(map(int, new[0:2]))))

            for iter_var, new in enumerate(new_temp):
                if new in old_temp:
                    to_dels.append(iter_var)

            for iter_var, to_del in enumerate(to_dels):
                fin_del = to_del - iter_var
                del new_data[fin_del]

            final_data = []
            final_temp = []

            for single_tors in new_data:
                if single_tors[0:2] not in final_temp:
                    final_data.append(single_tors)
                final_temp.append(single_tors[0:2])

            for final_tors in final_data:
                data.append(final_tors)

    return data


def file_write(final_file=None, final_data=None, case=None):

    """
    Takes final data to be written to xyz-File for Fortran conformer search, checks for the relevant keywords and adds
    all the $end statements.

    :param final_file: XYZ-File that is read in by Fortran conformer search.
    :param final_data: Data with relevant keywords and all ids added, redundancies and linear dependencies removed.
    :param case: Case to determine the type of reaction: intramol, bimol, shuttle.
    :return: Nothing is returned by this function, instead the final input-file for the conformer search is written.

    Input_gen and file_write so seemingly complicated to cope with all possible cases.
    Should be replaced in the future by much simpler direct argument passing to Fortran subroutines...
    """

    ind1 = 0
    ind2 = 0
    ind3 = 0
    ind4 = 0

    if [] in final_data:
        final_data.remove([])

    if ['$fadd'] not in final_data and case != 'intramol':
        final_data.insert(0, ['$fadd'])

    if ['$breaks'] not in final_data:
        final_data.append(['$breaks'])
        final_data.append([0, 0])

    for counter1, data_line in enumerate(final_data):
        if data_line[0] == '$fadd' and counter1 != 0:
            ind1 = counter1
    try:
        final_data.insert(ind1, ['$end'])
    except UnboundLocalError:
        pass

    for counter2, data_line in enumerate(final_data):
        if data_line[0] == '$breaks' and counter2 != 0:
            ind2 = counter2
    try:
        final_data.insert(ind2, ['$end'])
    except UnboundLocalError:
        pass

    for counter3, data_line in enumerate(final_data):
        if data_line[0] == '$tors' and counter3 != 0:
            ind3 = counter3
    try:
        final_data.insert(ind3, ['$end'])
    except UnboundLocalError:
        pass

    for counter4, data_line in enumerate(final_data):
        if data_line[0] == '$add' and counter4 != 0:
            ind4 = counter4
    try:
        final_data.insert(ind4, ['$end'])
    except UnboundLocalError:
        pass

    if final_data[-1] != ['$end']:
        final_data.append(['$end'])

    # Evaluation of number of conformers. Reduction one after another if too high. Starting from the bonds that are the
    # furthest apart from the reactive centers.

    tors_now = False
    rotations = []
    rot_ids = []
    fadd_rot = 0

    for id_count, list_item in enumerate(final_data):
        if list_item[0] == '$fadd':
            fadd_rot = final_data[id_count+1][2]
        if list_item[0] == '$tors':
            tors_now = True
        if tors_now:
            rotations.append(final_data[id_count])
            rot_ids.append(id_count)
        if list_item[0] == '$end':
            tors_now = False

    rotations = rotations[1:-1]
    rot_ids = rot_ids[1:-1]

    rot_factors = [360/rotation[2] for rotation in rotations]

    if fadd_rot > 0:
        rot_factors.insert(0, 360/fadd_rot)

    num_conf = np.prod(np.array(rot_factors))

    num_thresh = 5E8
    current_rots = []

    for rot_id in rot_ids:
        current_rots.append(final_data[rot_id][2])

    if num_conf > num_thresh:
        too_many = True

        while too_many:

            if num_conf < num_thresh:
                break

            elif all(current_rot == 120 or current_rot == 180 for current_rot in current_rots):

                id_shift = 0
                reduced = False
                rot_ids_temp = rot_ids[::-1].copy()

                for rot_id in rot_ids[::-1]:

                    if rot_id % 2 == 0:
                        final_data.pop(rot_id)
                        rot_ids_temp.remove(rot_id)
                        id_shift += 1

                    num_conf = 3**len(rot_ids_temp)

                    if num_conf < num_thresh:
                        reduced = True
                        break

                if reduced:
                    break

            else:
                current_rots = []

                for rot_id in rot_ids:
                    current_rots.append(final_data[rot_id][2])

                else:
                    for rot_id in rot_ids[::-1]:
                        if int(final_data[rot_id][2]) == 30:
                            final_data[rot_id][2] = final_data[rot_id][2]*2
                            num_conf /= 2
                            break

                        elif int(final_data[rot_id][2]) == 60:
                            final_data[rot_id][2] = final_data[rot_id][2]*2
                            num_conf /= 2
                            break

                        elif int(final_data[rot_id][2]) == 90:
                            final_data[rot_id][2] = final_data[rot_id][2]*2
                            num_conf /= 2
                            break

        with open(final_file, "a") as fileW:
            fileW.write('\n')
            first_join = ['\t'.join([str(x) if type(x) == int or type(x)==float else x for x in lst]) for lst in final_data]
            second_join = '\n'.join(first_join)
            fileW.write(second_join)
    else:
        with open(final_file, "a") as fileW:
            fileW.write('\n')
            first_join = ['\t'.join([str(x) if type(x) == int or type(x)==float else x for x in lst]) for lst in final_data]
            second_join = '\n'.join(first_join)
            fileW.write(second_join)

    return


def add_rule_stand(structures=None, rule=None):

    """
    Takes the ids for the add-rules and sorts them:
    - In general always first reactant first: [1, 0, 17, 3] -> [0, 1, 3, 17]
    - For two adds always intermolecular before intramolecular
    - ID-Shift for second reactant is removed

    :param structures: PMG-objects of the structures to be considered.
    :param rule: Reaction rule that will be standardized.
    :return: index_list: Single list (length 4 or 8) containing the shifted and appropriately sorted add-indices.
    """

    ruleAdd = rule['ADD']

    if len(structures) == 2:
        #id0 = len(structures[0])
        id0 = structures[0]["Number_of_Atoms"]
        index_list = [0, 1]

        ### One Add

        if len(ruleAdd) == 1:

            if ruleAdd[0][0] == 0 and ruleAdd[0][1] == 1:
                index_list.extend([ruleAdd[0][2], ruleAdd[0][3] - id0])

            if ruleAdd[0][0] == 1 and ruleAdd[0][1] == 0:
                index_list.extend([ruleAdd[0][3], ruleAdd[0][2] - id0])

        ### Two Adds

        elif len(ruleAdd) == 2:

            ### Both intermolecular

            if ruleAdd[0][0] == 0 and ruleAdd[0][1] == 1 and ruleAdd[1][0] == 0 and ruleAdd[1][1] == 1:
                index_list.extend([ruleAdd[0][2], ruleAdd[0][3] - id0, 0, 1, ruleAdd[1][2], ruleAdd[1][3] - id0])

            elif ruleAdd[0][0] == 0 and ruleAdd[0][1] == 1 and ruleAdd[1][0] == 1 and ruleAdd[1][1] == 0:
                index_list.extend([ruleAdd[0][2], ruleAdd[0][3] - id0, 0, 1, ruleAdd[1][3], ruleAdd[1][2] - id0])

            elif ruleAdd[0][0] == 1 and ruleAdd[0][1] == 0 and ruleAdd[1][0] == 0 and ruleAdd[1][1] == 1:
                index_list.extend([ruleAdd[0][3], ruleAdd[0][2] - id0, 0, 1, ruleAdd[1][2], ruleAdd[1][3] - id0])

            elif ruleAdd[0][0] == 1 and ruleAdd[0][1] == 0 and ruleAdd[1][0] == 1 and ruleAdd[1][1] == 0:
                index_list.extend([ruleAdd[0][3], ruleAdd[0][2] - id0, 0, 1, ruleAdd[1][3], ruleAdd[1][2] - id0])

            ### First one intermolecular

            elif ruleAdd[0][0] == 0 and ruleAdd[0][1] == 1 and ruleAdd[1][0] == 0 and ruleAdd[1][1] == 0:
                index_list.extend([ruleAdd[0][2], ruleAdd[0][3] - id0, 0, 0, ruleAdd[1][2], ruleAdd[1][3]])

            elif ruleAdd[0][0] == 0 and ruleAdd[0][1] == 1 and ruleAdd[1][0] == 1 and ruleAdd[1][1] == 1:
                index_list.extend(
                    [ruleAdd[0][2], ruleAdd[0][3] - id0, 1, 1, ruleAdd[1][2] - id0, ruleAdd[1][3] - id0])

            elif ruleAdd[0][0] == 1 and ruleAdd[0][1] == 0 and ruleAdd[1][0] == 0 and ruleAdd[1][1] == 0:
                index_list.extend([ruleAdd[0][3], ruleAdd[0][2] - id0, 0, 0, ruleAdd[1][2], ruleAdd[1][3]])

            elif ruleAdd[0][0] == 1 and ruleAdd[0][1] == 0 and ruleAdd[1][0] == 1 and ruleAdd[1][1] == 1:
                index_list.extend(
                    [ruleAdd[0][3], ruleAdd[0][2] - id0, 1, 1, ruleAdd[1][3] - id0, ruleAdd[1][2] - id0])

            ### Second one intermolecular

            elif ruleAdd[0][0] == 0 and ruleAdd[0][1] == 0 and ruleAdd[1][0] == 0 and ruleAdd[1][1] == 1:
                index_list.extend([ruleAdd[1][2], ruleAdd[1][3] - id0, 0, 0, ruleAdd[0][2], ruleAdd[0][3]])

            elif ruleAdd[0][0] == 0 and ruleAdd[0][1] == 0 and ruleAdd[1][0] == 1 and ruleAdd[1][1] == 0:
                index_list.extend([ruleAdd[1][3], ruleAdd[1][2] - id0, 0, 0, ruleAdd[0][2], ruleAdd[0][3]])

            elif ruleAdd[0][0] == 1 and ruleAdd[0][1] == 1 and ruleAdd[1][0] == 0 and ruleAdd[1][1] == 1:
                index_list.extend(
                    [ruleAdd[1][2], ruleAdd[1][3] - id0, 1, 1, ruleAdd[0][2] - id0, ruleAdd[0][3] - id0])

            elif ruleAdd[0][0] == 1 and ruleAdd[0][1] == 1 and ruleAdd[1][0] == 1 and ruleAdd[1][1] == 0:
                index_list.extend(
                    [ruleAdd[1][3], ruleAdd[1][2] - id0, 1, 1, ruleAdd[0][2] - id0, ruleAdd[0][3] - id0])

            else:
                print('Error in rule-Handling')
                pass

    else:
        for i in range(len(ruleAdd)):
            ruleAdd[i][2:4] = sorted(ruleAdd[i][2:4])
        index_list = [item for sublist in ruleAdd for item in sublist]

    return index_list


def break_rule_stand(structures=None, rule=None):

    """
    Breaks not as crucial as adds, therefore function not as sophisticated. If appropriate then index-shift is removed.

    :param structures: PMG-objects of the structures to be considered.
    :param rule: Break rule that will be standardized.
    :return: index_list: Single list (length 3 or 6) containing the shifted break-indices. If no breaks then None is
    returned.
    """

    ruleBreak = rule['BREAK']

    if len(ruleBreak) > 0:
        index_list = []
        for i in range(len(ruleBreak)):
            ruleBreak[i][1] = ruleBreak[i][1] - ruleBreak[i][0] * len(structures[0])
            ruleBreak[i][2] = ruleBreak[i][2] - ruleBreak[i][0] * len(structures[0])
            index_list.extend(ruleBreak[i])
    else:
        index_list = []

    return index_list


def genDisConn(structures=None):

    """
    This function generates a connectivity table for the disconnected reactants and writes it into the file
    'dis_conn.dat'. The fragmented structure is also written out as 'temp_precomp.xyz'.

    :param structures: PMG-objects of the structures to be considered.
    :return: None.
    """
    if len(structures) > 1:
        temp = structures[0].copy()
        if len(structures) == 1:
            pass
        elif len(structures) == 2:
            temp2 = structures[1].copy()
            temp2.translate_sites(vector=[0, 0, 1000])

            for i in temp2:
                temp.append(species=i.species_string, coords=i.coords, validate_proximity=False)

        temp.to(fmt='xyz', filename='temp_precomp.xyz')

    return


def reg_precomp(preCompList=None, input_name=''):

    """
    Function to properly register the generated pre-complex. If the conformer search/geometry optimization was
    successful then the conf*.xyz - File is read in, if not, the initially aligned structure. One can also pass a
    specific filename to the function.

    :param preCompList: List of pre-complexes where the generated structure will be appended to upon registration.
    :param input_name: Optional Filename-argument (without .xyz-extension).
    :return: None.
    """

    conf_create = False
    for conf_out in os.listdir('.'):
        if fnmatch.fnmatch(conf_out, 'conf0*.xyz'):
            found_conf = Molecule.from_file(conf_out)
            conf_create = True
        elif fnmatch.fnmatch(conf_out, 'precomp_pre_conf.xyz'):
            pre_conf = Molecule.from_file(conf_out)
        elif fnmatch.fnmatch(conf_out, input_name + '.xyz'):
            pre_conf = Molecule.from_file(conf_out)

    if conf_create:
        preCompList.append(found_conf)

    elif not conf_create:
        preCompList.append(pre_conf)

    return


def project_onto_plane(v=None, n=None):

    """
    Function that projects a given vector v onto a plane defined by its perpendicular normal vector.

    :param x: Vector that should be projected onto the plane.
    :param n: Normal vector of the plane
    :return:
    """

    d = np.dot(v, n) / np.linalg.norm(n)
    p = [d * norm_vec(n)[i] for i in range(len(n))]
    proj_vec = np.asarray([v[i] - p[i] for i in range(len(v))])

    return proj_vec


def dihed_avg(atomCoord=None, neighCoord=None):

    """
    Function that calculates the average dihedral angle of the four neighbouring atoms, or the angle between the
    neighbouring atoms and the reactive atom, respectively. This leads to a much more reliable classification in the
    topology analysis instead of just calculating one random angle.

    In principal one should also implement this in Fortran where so far only one angle is considered. Due to the
    higher effort that would have to be undertaken as a functionality like 'itertools.permutations' is not available
    there, this is only done in Python so far. For the alignment of the second add vectors in the Fortran code the
    choice of a larger angle range is not as problematic and the importance of the topology analysis in general also is
    much lower.

    :param atomCoord: Coordinates of reactive atom.
    :param neighCoord: Coordinates of neighbouring atoms (3 or 4).
    :return: dih_avg: Average dihedral angle between the four neighbouring atoms or the three neighbouring and the
                      reactive atom.
    """

    dihed_permuts = None

    dihed_sum = 0
    if len(neighCoord) == 4:
        dihed_permuts = list(itertools.permutations([0, 1, 2, 3]))
        for dihed_permut in dihed_permuts:
            dihed_sum += np.absolute(dihedral_func(np.array([neighCoord[dihed_permut[i]] for i in range(4)])))

    elif len(neighCoord) == 3:
        dihed_permuts = list(itertools.permutations([0, 1, 2]))
        for dihed_permut in dihed_permuts:
            dihed_sum += np.absolute(dihedral_func(np.array([neighCoord[dihed_permut[0]],
                                                             atomCoord,
                                                             neighCoord[dihed_permut[1]],
                                                             neighCoord[dihed_permut[2]]])))

    if dihed_permuts != None:
        dihed_avg = np.absolute(dihed_sum / len(dihed_permuts))

    return dihed_avg


def angle_avg(atomCoord=None, neighCoord=None):

    """
    Function that calculates the average angle between the reactive and its neighbouring atoms. This leads to a much
    more reliable classification in the topology analysis instead of just calculating one random angle.

    In principal one should also implement this in Fortran where so far only one angle is considered. Due to the
    higher effort that would have to be undertaken as a functionality like 'itertools.permutations' is not available
    there, this is only done in Python so far. For the alignment of the second add vectors in the Fortran code the
    choice of a larger angle range is not as problematic and the importance of the topology analysis in general also is
    much lower.

    :param atomCoord: Coordinates of reactive atom.
    :param neighCoord: Coordinates of neighbouring atoms (2 to 4).
    :return: ang_avg: Average angle between the reactive and its neighbouring atoms.
    """

    angle_permuts = None
    ang_avg = None

    angle_sum = 0
    if len(neighCoord) == 4:
        angle_permuts = list(itertools.permutations([0, 1, 2, 3]))

    elif len(neighCoord) == 3:
        angle_permuts = list(itertools.permutations([0, 1, 2]))

    elif len(neighCoord) == 2:
        angle_permuts = list(itertools.permutations([0, 1]))

    if angle_permuts != None:
        for angle_permut in angle_permuts:
            angle_sum += np.absolute(angle_func((atomCoord - neighCoord[angle_permut[0]]),
                                                (atomCoord - neighCoord[angle_permut[1]])))

        ang_avg = np.absolute(angle_sum / len(angle_permuts))

    return ang_avg


def get_key(inp_dict=None, ra=None):

    """
    Small helper function that retrieves the relevant topology dictionary key from a reactive atom ID.
    It is important to mention that the IDs in the Dictionary and of the reactive atom here are the Fortran-ones and
    therefore start at 1, NOT at 0.

    :param inp_dict: Topology dictionary set up as: str(ID + NN + type_abbr + atomic_sym): [rvecs]
    :param ra: Reactive atom.
    :return: Returns the 8 or 9-character string that functions as key for the topology dictionary to retrieve the
             rection vectors from it.
    """

    out_key = None
    for i in inp_dict:
        if int(i[:3]) == ra:
            out_key = i

    if out_key is None:
        print('Problem with retrieval of topology key from reactive atom ID.', ra)
        out_key = '0000xxxX'

    return out_key


def gen_sadd(path=None, sadd_factor=None):

    """
    Depending on the path length, define the threshold for the second add distance
    """

    if len(path) >= 2:
        sadd_thresh = sadd_factor * 4 * (1 + (1 / len(path)))
    else:
        sadd_thresh = 6

    return str(round(sadd_thresh, 2))

def readMols(name, type="xyz", workdir="."):

    """
    Read XYZ files and generate openbabel and pymatgen molecule objects.
    """

    datei = ".".join([name, type])
    try:
        pmg = Molecule.from_file(os.path.join(workdir, datei))
        obconversion = ob.OBConversion()
        obconversion.SetInFormat("xyz")
        mol = ob.OBMol()
        obconversion.ReadFile(mol, str(datei))
        obm = mol
    except Exception:
        try:
            exc_info = sys.exc_info()
            try:
                raise TypeError()
            except:
                pass
        finally:
            traceback.print_exception(*exc_info)
            del exc_info
    return (pmg, obm)

def check_for_rings(educt=None, at_num_1=None, at_num_2=None):

    found = False
    neighbors = False

    edges = connectivity(mol=educt)

    if len(edges) != 0:  # otherwise just an atom!

        sorted_list = sorted([at_num_1, at_num_2])
        if (sorted_list[0], sorted_list[1]) in edges:
            neighbors = True

        G = nx.Graph()
        G.add_edges_from(edges)
        rings = nx.cycle_basis(G)
        for ring in rings:
            if at_num_1 in ring and at_num_2 in ring:
                found = True

    return found, neighbors

def connectivity(mol=None):

    """
    Determine connectivity based on distance-cutoffs and covalent radii determination.
    """

    mol_ob = BabelMolAdaptor(mol).openbabel_mol
    edges = []
    symbols = []

    for atom in ob.OBMolAtomIter(mol_ob):
        # symbols.append(ob.OBElements.GetSymbol(atom.GetAtomicNum()))
        symbols.append(element_dict[atom.GetAtomicNum()])

    #Generate connectivity table:

    for i, atom1 in enumerate(ob.OBMolAtomIter(mol_ob)):
        for j, atom2 in enumerate(ob.OBMolAtomIter(mol_ob)):
            if i != j:
                symbol1 = symbols[i]
                symbol2 = symbols[j]
                coord1 = np.array([atom1.GetX(), atom1.GetY(), atom1.GetZ()])
                coord2 = np.array([atom2.GetX(), atom2.GetY(), atom2.GetZ()])
                dist = np.linalg.norm(coord1-coord2)
                radius = kov_rad(symbol1, symbol2)
                if (dist < radius and dist > 0.0):
                    if atom1.GetIdx() > atom2.GetIdx() and (atom2.GetIdx(), atom1.GetIdx()) not in edges:
                        edges.append((atom2.GetIdx(), atom1.GetIdx()))
                    elif atom2.GetIdx() > atom1.GetIdx() and (atom1.GetIdx(), atom2.GetIdx()) not in edges:
                        edges.append((atom1.GetIdx(), atom2.GetIdx()))

    return edges

def repair_OBmol_conn(obm=None):

    """
    Determine differences between the distance-based connectivity and the
    OB "interpreted" connectiviy.
    Modify an openbabel molecule object according to the distance-based connectivity.
    """
    edgesOB = []

    for bond in ob.OBMolBondIter(obm):
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetIdx() > atom2.GetIdx() and (atom2.GetIdx(), atom1.GetIdx()) not in edgesOB:
            edgesOB.append((atom2.GetIdx(), atom1.GetIdx()))
        elif atom2.GetIdx() > atom1.GetIdx() and (atom1.GetIdx(), atom2.GetIdx()) not in edgesOB:
            edgesOB.append((atom1.GetIdx(), atom2.GetIdx()))

    # get connectivitities based on distances
    edges_dist = sorted(connectivity(obm))

    # check difference to OBbond-connectivities
    diff_dist_OB = set(edges_dist) - set(edgesOB) # we have to add these connectivities to OBmol

    # modify OBmol object
    modifyOBM(obm=obm, diff=diff_dist_OB)


def modifyOBM(obm, diff):

    """
    Modify an openbabel molecule object according to the distance-based connectivity.
    """

    for item in diff:
        obm.AddBond(item[0], item[1], 1)

    return obm

def getEducts(mols, verbose=False, adds=None):
    """
    Get list of "educt dictionaries" containing a "name", "type", "obm" (openbabel molecule object),
    "pmg" (pymatgen molecule object) and "Number_of_Atoms"
    """
    educts = []
    for c, mol in enumerate(mols):
        if isinstance(mol["name"], ob.OBMol):
            educt = {"name": "m0" + str(c + 1), "type": "xyz"}
            obm = mol["name"]
            pmg = BabelMolAdaptor(obm).pymatgen_mol
        elif isinstance(mol["name"], Molecule):  # PyMatGen
            educt = {"name": "m0" + str(c + 1), "type": "xyz"}
            pmg = mol["name"]
            obm = BabelMolAdaptor(pmg).openbabel_mol
        else:
            name = mol["name"]
            type = mol["type"]
            educt = {"name": name, "type": type}
            pmg, obm = readMols(name)
        educt["Number_of_Atoms"] = len(pmg)
        educt["pmg"] = pmg
        educt["obm"] = obm

        # evaluate disconnections (e.g., Lewis-Acids)
        edges = []

        repair_OBmol_conn(educt["obm"])

        for bond in ob.OBMolBondIter(educt["obm"]):
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            edges.append((atom1.GetIdx(), atom2.GetIdx()))

        if len(edges) != 0: #otherwise just an atom!
            G = nx.Graph()
            G.add_edges_from(edges)
            if not nx.is_connected(G):
                if adds != None:
                    for add in adds:
                        if int(add[0]) == c:
                            educt["obm"].AddBond(add[1], add[2], 1)
                else:
                    sub_graphs = nx.connected_component_subgraphs(G)
                    graphs = [sg.nodes() for sg in sub_graphs]
                    dist_min = 100
                    atom1_min = 0
                    atom2_min = 0
                    for node1 in graphs[0]:
                        for node2 in graphs[1]:
                            vec1 = np.array([educt["obm"].GetAtom(node1).GetX(), educt["obm"].GetAtom(node1).GetY(),
                                             educt["obm"].GetAtom(node1).GetZ()])
                            vec2 = np.array([educt["obm"].GetAtom(node2).GetX(), educt["obm"].GetAtom(node2).GetY(),
                                             educt["obm"].GetAtom(node2).GetZ()])
                            dist = np.linalg.norm(vec1-vec2)
                            if dist < dist_min:
                                dist_min = dist
                                atom1_min = node1
                                atom2_min = node2
                    try:
                        add1 = atom1_min
                        add2 = atom2_min
                        educt["obm"].AddBond(add1, add2, 1)
                        if verbose:
                            print("##############################")
                            print("NOTE: ADDED bond:", add1, add2)
                            print("##############################")
                    except UnboundLocalError:
                        print("Addition of new bond was not possible!")

        educts.append(educt)

    return educts


def genDisConn(structures):

    """
    This function generates a connectivity table for the disconnected reactants and writes it into the file
    'dis_conn.dat'. The fragmented structure is also written out as 'temp_precomp.xyz'.

    :param structures: PMG-objects of the structures to be considered.
    :return: None.
    """
    if len(structures) > 1:
        temp = structures[0]["pmg"].copy()
        if len(structures) == 1:
            pass
        elif len(structures) == 2:
            temp2 = structures[1]["pmg"].copy()
            temp2.translate_sites(vector=[0, 0, 1000])

            for i in temp2:
                temp.append(species=i.species_string, coords=i.coords, validate_proximity=False)

        temp.to(fmt='xyz', filename='temp_precomp.xyz')

    return
