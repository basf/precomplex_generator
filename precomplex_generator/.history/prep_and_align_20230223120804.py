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

import numpy as np
import pymatgen as pmg
import openbabel as ob
from pymatgen.core.structure import Molecule
from get_all_torsions_mod import get_tors
from tools import file_read_dict, file_read_list, angle_func, norm_vec, center_of_points, rot_mat, project_onto_plane, dihed_avg, angle_avg, kov_radius, get_key
import os
import sys
import math
from copy import deepcopy
from element_dict import element_dict


def topo_analysis(educt=None, educt_name=None, own_angle=None, shuttle=False):

    """
    Function that carries out the topology analysis of the reactants. First, the connectivity table will be created and
    the rotatable bonds identified by the Fortran-Routines. It is then iterated through the atoms and depending on the
    number of neighbours the relevant dihedral and neighbour-angles are calculated. These information then define the
    geometric case which the atom belongs to and how the reaction vector is defined. These are then written into a dic-
    tionary where the key is composed of:

    - A three-character string representing the ID from the Fortran code (starting at 1 instead of 0), e.g. '004'.
    - The number of neighbours as a one-character string.
    - A three-character abbreviation for the geometry case the atom belongs to.
    - And the chemical symbol which can be either one or two characters long and is therefore put at the end to easily
      access it by the string slice [7:]

    If there is more than one possible definition for the angles, the average over all combinations of neighbours is
    calculated. This allows also skewed structures to be correctly identified, removes any bias by the ordering of the
    xyz-File and lowers the needed range significantly. If the atom cannot be assigned to any of the cases a warning is
    printed out, the reaction vector is defined along the z-axis and the type abbreviation is set to 'xxx'.

    :param educt_name: names of reactant structure files (without .xyz ending)
    :param educt: dictionary containing pymatgen and openbabel molecule objects of reactant structures
    :return: dictionary with all collected information ("analysis")
    """

    analysis = {}
    angle = 0
    dihed = 0
    short_type = 'xxx'

    if educt_name == 'precomp_pre_solv':

        xyz_list = file_read_list(filename='precomp_pre_solv.xyz')
        conn_dict = file_read_dict(filename='../../../dis_conn.dat')

    else:
        write_out_conn(educt)
        get_tors(educt_name + '_conn.dat')
        os.rename('tors.dat', educt_name + '_tors.dat')

        analysis = {}
        xyz_list = file_read_list(filename=educt_name + '.xyz')
        conn_dict = file_read_dict(educt_name + '_conn.dat')

    x_sum, y_sum, z_sum = 0, 0, 0

    for xyz_coord in xyz_list:
        if xyz_coord == []:
            continue
        if xyz_coord[0][0] == "$":
            break
        x_sum += float(xyz_coord[1])
        y_sum += float(xyz_coord[2])
        z_sum += float(xyz_coord[3])

    cop_all = np.asarray([x_sum, y_sum, z_sum])/len(xyz_list)
    educt_pmg = Molecule.from_file(educt_name+'.xyz')
    com_all = educt_pmg.center_of_mass

    copm = 0.5*cop_all + 0.5*com_all

    for conn_key in conn_dict.keys():
        numNeigh = len(conn_dict[conn_key])
        symbol = xyz_list[int(conn_key) - 1][0]
        atomCoord = np.asarray(list(map(float, xyz_list[int(conn_key) - 1][1:])))
        neighCoord = []
        neighCoord_new = [deepcopy(atomCoord)]

        for j in conn_dict[conn_key]:
            neighCoord_new.append(np.asarray(
                list(map(float, xyz_list[int(j) - 1][1:]))))
            neighCoord.append(np.asarray(
                list(map(float, xyz_list[int(j) - 1][1:]))))

        cop_n = center_of_points(array=np.asarray(neighCoord_new))

        if numNeigh > 4:
            sys.exit("Error: It is not yet possible to define a reaction vector for an atom with more than 4 neighbors!")

        elif numNeigh == 4:

            angle = angle_avg(atomCoord=atomCoord, neighCoord=neighCoord)
            dihed = dihed_avg(atomCoord=atomCoord, neighCoord=neighCoord)

            if 100 <= angle <= 120 and 60 <= dihed <= 80:
                rvecs = [norm_vec(atomCoord - neighCoord[0]), norm_vec(atomCoord - neighCoord[1]),
                         norm_vec(atomCoord - neighCoord[2]), norm_vec(atomCoord - neighCoord[3])]

                copm_angles = {}
                for neigh_i in range(len(neighCoord)):
                    copm_angles[neigh_i] = np.absolute(angle_func(norm_vec(atomCoord - neighCoord[neigh_i]), copm))

                short_type = 'tet'

                rvecs_small_ring = rvecs

        elif numNeigh == 3:

            angle = angle_avg(atomCoord=atomCoord, neighCoord=neighCoord)
            dihed = dihed_avg(atomCoord=atomCoord, neighCoord=neighCoord)

            if 110 <= angle <= 130 and (dihed <= 10 or 160 <= dihed <= 190):
                v1 = atomCoord - neighCoord[0]
                v2 = atomCoord - neighCoord[1]
                rvecs = [norm_vec(np.cross(v1, v2)), norm_vec(np.cross(v2, v1))]
                short_type = 'tpl'

            elif 80 <= angle <= 130 and 80 <= dihed <= 160:
                rvecs = [norm_vec(atomCoord - cop_n)]
                short_type = 'tpy'

            rvecs_small_ring = rvecs

        elif numNeigh == 2:
            angle = angle_avg(atomCoord=atomCoord, neighCoord=neighCoord)
            if 150 <= angle <= 190:
                v_plane = neighCoord[0] - neighCoord[1]  # Orthogonal vector defining plane (normal vector)
                rvecs = [-norm_vec(project_onto_plane((copm-atomCoord), v_plane))]
                short_type = 'lin'

            elif angle < 150:
                v1 = atomCoord - neighCoord[0]
                v2 = atomCoord - neighCoord[1]
                rvecs_horizontal = [norm_vec(np.cross(v1, v2)), norm_vec(np.cross(v2, v1))]
                rvecs_start = [norm_vec(atomCoord - cop_n)]
                rvecs = [norm_vec(rvecs_start[0] + rvecs_horizontal[0])]
                short_type = 'ang'

            rvecs_small_ring = rvecs

        elif numNeigh == 1:

            if own_angle != None and (symbol != "H" or shuttle):
                own_ang = (180 - own_angle) * math.pi / 180
                v1 = atomCoord - neighCoord[0]
                v2 = copm - neighCoord[0]
                vec = norm_vec(atomCoord - neighCoord[0])

                # Comment: if atom coordinates are 0 0 0, we get a problem for the axis, when we divide by 0!
                # ToDo: shifting by 0 0 0.01 helps, but then a 180 angle and not the own_angle is used!

                cp = np.cross(v1, v2) + np.array([0, 0, 1E-2])

                """Rodrigues rotation formula"""
                axis = cp / np.linalg.norm(cp)
                term1 = vec * np.cos(own_ang)
                term2 = (np.cross(axis, vec)) * np.sin(own_ang)
                term3 = axis * ((1 - np.cos(own_ang)) * axis.dot(vec))
                rvecs = [ norm_vec(term1 + term2 + term3) ]
                short_type = 'end'
            else:
                rvecs = [norm_vec(atomCoord - neighCoord[0])]
                short_type = 'end'

            # for small rings, assume 180°, even if own_angle is different!
            rvecs_small_ring = [norm_vec(atomCoord - neighCoord[0])]

        elif numNeigh == 0:
            rvecs = [np.array([0, 0, 1])]
            rvecs_small_ring = rvecs
            short_type = 'ato'

        try:
            analysis[str(conn_key).zfill(3) + str(numNeigh) + short_type + symbol] = {"normal": rvecs, "small_rings": rvecs_small_ring}
        except UnboundLocalError:
            print('\nAssignment of atom {0} to a geometry case failed! The values are:'.format(
                   str(conn_key).zfill(3) + symbol))
            print('\n{0} (number of neighbours), {1} (angle), {2} (dihedral angle).'.format(
                   numNeigh, angle, dihed))
            print('\nGeometry might be quite skewed and range in first_builder.py should be increased!')
            print('The unit vector along the z-axis will be set as reaction vector for this atom.\n')
            analysis[str(conn_key).zfill(3) + str(numNeigh) + short_type + symbol] = None

    return analysis


def write_out_conn(educt=None):

    """
    Write out connectivity file.
    """

    name = educt["name"]+"_conn.dat"
    neighbour_dic = {}

    for atom in ob.OBMolAtomIter(educt["obm"]):
        id1 = atom.GetIdx()
        # symbol1 = ob.OBElementTable().GetSymbol(atom.GetAtomicNum())
        symbol1 = element_dict[atom.GetAtomicNum()]
        tmp = []
        nneigh = 0
        for neigh in ob.OBAtomAtomIter(atom):
            nneigh += 1
            id2 = neigh.GetIdx()
            neighbour_string = str(id2)
            tmp.append(neighbour_string)
        neighbour_dic[str(id1)] = str(id1) + " " + symbol1 + " " + str(nneigh)+" "+" ".join(sorted(tmp))

    with open(name, "w") as out:
        for key in neighbour_dic.keys():
            out.write(neighbour_dic[key]+"\n")


def write_out_disconn(educts=None):

    """
    Write out connectivity of disconnected molecules (in bimolecular case).
    """

    name = "dis_conn.dat"
    if os.path.isfile(name):
        os.remove(name)

    neighbour_dic = {}
    shift = educts[0]["Number_of_Atoms"]

    for i, mol in enumerate(educts):
        indexShift = i*shift
        for atom in ob.OBMolAtomIter(mol["obm"]):
            id1 = atom.GetIdx()+indexShift
            # symbol1 = ob.OBElementTable().GetSymbol(atom.GetAtomicNum())
            symbol1 = element_dict[atom.GetAtomicNum()]
            tmp = []
            nneigh = 0
            for neigh in ob.OBAtomAtomIter(atom):
                nneigh += 1
                id2 = neigh.GetIdx()+indexShift
                neighbour_string = str(id2)
                tmp.append(neighbour_string)
            neighbour_dic[str(id1)] = str(id1) + " " + symbol1 + " " + str(nneigh) + " " + " ".join(sorted(tmp))

    with open(name, "w") as out:
        for key in neighbour_dic.keys():
            out.write(neighbour_dic[key] + "\n")

def addSolv(pmgs, solv_acc, cplx_don, topos, shuttle_tup, prev_ras, len_first, len_both, intelligentDistance):

    prev_ra1 = prev_ras[2] + len_first * prev_ras[0]
    prev_ra2 = prev_ras[3] + len_first * prev_ras[1]

    cplx, solv = pmgs[0], pmgs[1]

    cen_cplx = pmg.core.operations.SymmOp.from_rotation_and_translation(translation_vec=-cplx[cplx_don-1].coords)
    cplx.apply_operation(cen_cplx)
    cen_solv = pmg.core.operations.SymmOp.from_rotation_and_translation(translation_vec=-solv[solv_acc-len_both-1].coords)
    solv.apply_operation(cen_solv)

    topos['cplx'] = topo_analysis(educt_name='precomp_pre_solv')

    cplx_keys = topos['cplx'].keys()
    for cplx_key in cplx_keys:
        if int(cplx_key[0:3]) == cplx_don:
            rvec_cplx = np.asarray(topos['cplx'][cplx_key]["normal"][0])
            rad_cplx = kov_radius(cplx[int(cplx_key[0:3]) - 1].species_string)

    solv_keys = topos[shuttle_tup[0]]
    for solv_key in solv_keys:
        if (solv_acc-len_both) == int(solv_key[0:3]):
            rvec_solv = np.asarray(topos[shuttle_tup[0]][solv_key]["normal"][0])
            rad_solv = kov_radius(solv[int(solv_key[0:3]) - 1].species_string)

    solv_copy, cplx_copy = solv.copy(), cplx.copy()

    R_align = rot_mat(-rvec_solv, rvec_cplx)

    z = pmg.core.operations.SymmOp.from_rotation_and_translation(rotation_matrix=R_align)
    solv_copy.apply_operation(z)

    # First translation of solvent into direction of reaction vector of complex
    t2 = pmg.core.operations.SymmOp.from_rotation_and_translation(
        translation_vec = 0.5*intelligentDistance * (rad_solv + rad_cplx) * rvec_cplx + np.array([0, 0, 1E-3]))
    solv_copy.apply_operation(t2)

    # Calculation of tilt-vector and matrix
    fa_vec = -norm_vec(cplx_copy[prev_ra2 - 1].coords - cplx_copy[prev_ra1 - 1].coords)
    hbond_vec = norm_vec(solv_copy[solv_acc - len_both - 1].coords - cplx_copy[cplx_don - 1].coords)

    #tilt_vec = norm_vec(fa_vec + 3*hbond_vec) # JG original
    tilt_vec = norm_vec(fa_vec + 1 * hbond_vec) # the smaller the prefactor, the larger the tilt!!!

    R_tilt = rot_mat(rvec_cplx, tilt_vec)

    t3 = pmg.core.operations.SymmOp.from_rotation_and_translation(
        translation_vec = 0.5 * (rad_solv + rad_cplx) * tilt_vec + np.array([0, 0, 1E-3]), rotation_matrix=R_tilt)
        # translation_vec = 0.5 * (rad_solv + rad_cplx) * norm_vec(rvec_cplx) + np.array([0, 0, 1E-3]), rotation_matrix=R_tilt)
    solv_copy.apply_operation(t3)

    for i in solv_copy:
        cplx_copy.append(species=i.species_string, coords=i.coords, validate_proximity=False)

    cplx_copy.to(filename='cplx_with_solv.xyz')

    return cplx_copy

def rvec_align(educt_names=None, educts=None, addRule=None, breakRules=None, topos=None, user_dist=None):

    """
    Funtion that aligns the reactants based on a reaction rule and the reaction vectors of the add-atoms. The names and
    pmg-Molecule objects of the educts, the adds and breaks for the relevant rule, the topology dictionary and the
    alignment-distance factor are needed.

    :param educt_names: String containing the names of the educts.
    :param pmgs: Pymatgen Molecule objects of the educts.
    :param addRule: Standardized Adds of the reaction rule.
    :param breakRules: Standardized Breaks of the reaction rules.
    :param topos: Dictionary containing the topology dictionaries for the educts with their names as keys.
    :param distance: Distance factor that can be changed by the user (default: 1).
    :return: Returns a list of pymatgen objects with all the geometric possibilities (1,2 or 4) of the pre-complexes.

    The largest and most complicated part of the code is needed to identify the relevant break-rules for tetrahedral
    atoms that are involved in bond-formations. If this is not the case, the routine is quite simple.

    """

    # Setting up relevant variables

    ra1 = addRule[2]
    ra2 = addRule[3]

    mol1 = educts[0]["pmg"].copy()
    mol2 = educts[1]["pmg"].copy()

    topo1_keys = topos[educt_names[0]].keys()
    topo2_keys = topos[educt_names[1]].keys()
    key1 = get_key(inp_dict=topo1_keys, ra=ra1)
    key2 = get_key(inp_dict=topo2_keys, ra=ra2)
    topo_keys = [key1, key2]

    problem = False

    if topos[educt_names[0]][key1] != None:
        rvecs1 = topos[educt_names[0]][key1]["normal"]
    else:
        problem = True # assignment of reaction vector did not work!

    if topos[educt_names[1]][key2] != None:
        rvecs2 = topos[educt_names[1]][key2]["normal"]
    else:
        problem = True

    rad1 = kov_radius(str(mol1[ra1 - 1].species_string))
    rad2 = kov_radius(str(mol2[ra2 - 1].species_string))
    sum_rads = rad1+rad2
    precplx = []

    # Centering both reactants on their reactive atoms. This is crucial for the following rotation and translation
    cen1 = pmg.core.operations.SymmOp.from_rotation_and_translation(translation_vec=-mol1[ra1-1].coords)
    mol1.apply_operation(cen1)

    cen2 = pmg.core.operations.SymmOp.from_rotation_and_translation(translation_vec=-mol2[ra2-1].coords)
    mol2.apply_operation(cen2)

    used_rvecs = []

    if problem:
        print('Error in RVEC-alignment function...')
        precplx.append(None)

    # Alignment of reaction vectors if no tetrahedral species are involved
    elif (len(rvecs1) < 4 and len(rvecs2) < 4) or ((((len(rvecs1) == 4 and key1[7:] != 'C') or
                                                   (len(rvecs2) == 4) and key2[7:] != 'C')) and len(breakRules) == 0):

        align_clash = False

        for rvec1 in rvecs1:
            for rvec2 in rvecs2:

                used_rvecs_dict = {0:{ra1:[]}, 1:{ra2:[]}} # the already used reaction vectors from educt 0 and 1

                mol1temp = mol1.copy()

                R = rot_mat(-rvec1, rvec2)

                z = pmg.core.operations.SymmOp.from_rotation_and_translation(rotation_matrix=R)
                mol1temp.apply_operation(z)

                t2 = pmg.core.operations.SymmOp.from_rotation_and_translation(
                    translation_vec=user_dist * sum_rads * rvec2 + np.array([0, 0, 1E-2]))

                mol1temp.apply_operation(t2)

                moltemp_test = mol1temp.copy()
                moltemp_far = mol1temp.copy()

                for i in mol2:
                    mol1temp.append(species=i.species_string, coords=i.coords, validate_proximity=False)

                    try:
                        moltemp_test.append(species=i.species_string, coords=i.coords, validate_proximity=True)
                    except ValueError:
                        align_clash = True

                precplx.append(mol1temp)

                used_rvecs_dict[0][ra1] = rvec1
                used_rvecs_dict[1][ra2] = rvec2

                used_rvecs.append(used_rvecs_dict)

                if align_clash:

                    t_far = pmg.core.operations.SymmOp.from_rotation_and_translation(
                        translation_vec = 10 * sum_rads * np.array([0, 0, 1]))
                    moltemp_far.apply_operation(t_far)

                    for i in mol2:
                        moltemp_far.append(species=i.species_string, coords=i.coords, validate_proximity=False)

                    moltemp_far.to(filename="far_sep.xyz")


    elif (len(rvecs1) == 4 or len(rvecs2) == 4) and len(breakRules) > 0:


        pmg_dict = {1: mol1, 2: mol2}
        ra_dict = {1: ra1, 2: ra2}
        rvec_dict = {1: rvecs1, 2: rvecs2}

        sn2_struc = sn2_align(pmg_dict=pmg_dict, ra_dict=ra_dict, sum_rads=sum_rads, breakRules=breakRules,
                              addRule=addRule, rvec_dict=rvec_dict, user_dist=user_dist)

        precplx.append(sn2_struc)

    elif (len(rvecs1) == 4 or len(rvecs2) == 4) and len(breakRules) == 0:
        print("NOT yet implemented: addition to tetrahedral without break!")
        sys.exit()

    else:
        print('Error in RVEC-alignment function...')
        precplx.append(None)

    return precplx, used_rvecs, topo_keys


def sn2_align(pmg_dict, ra_dict=None, sum_rads=None, addRule=None, breakRules=None, rvec_dict=None, user_dist=None):

    """
    Take care of SN2-like reactions (= addition to tetrahedral center with simultaneous breakage of one bond
    connected to the "center" atom)
    :param pmg_dict:
    :param ra_dict:
    :param sum_rads:
    :param addRule:
    :param breakRules:
    :param rvec_dict:
    :param user_dist:
    :return: "return_mol"= ...
    """

    num_r1 = len(pmg_dict[1])

    ra1_tot = ra_dict[1] + addRule[0] * len(pmg_dict[1])
    ra2_tot = ra_dict[2] + addRule[1] * len(pmg_dict[1])
    ra_tot_dict = {1: ra1_tot, 2: ra2_tot} # reactive atoms in "summed up" atom-number format

    direction = 1
    invert = False

    if len(breakRules) == 3:  # ONE break rule

        cl_break = [breakRules[0]]

        for i in range(1, 3):
            breakRules[i] = breakRules[i] + breakRules[0] * num_r1

        if len(rvec_dict[1]) == 4 and len(rvec_dict[2]) < 4: # if len(rvec_dict) ==4 --> this is a tetrahedral env.
            tet = 1; oth = 2

        elif len(rvec_dict[1]) < 4 and len(rvec_dict[2]) == 4:
            tet = 2; oth = 1
            invert = True

        for x in breakRules[1:3]: # loop over atoms in break rule (summed up)
            if x is not ra_tot_dict[tet]:
                cl_break.append(x)
        if len(cl_break) > 2:
            print("Add and break are NOT connected!!!")

    elif len(breakRules) == 6: # TWO break rules

        first_br = breakRules[0:3]
        second_br = breakRules[3:6]

        for i in range(1, 3):
            first_br[i] = first_br[i] + first_br[0] * num_r1
            second_br[i] = second_br[i] + second_br[0] * num_r1

        if len(rvec_dict[1]) == 4 and len(rvec_dict[2]) < 4:
            tet = 1
            oth = 2

        elif len(rvec_dict[1]) < 4 and len(rvec_dict[2]) == 4:
            tet = 2
            oth = 1
            invert = True

        else:
            tet = 1
            oth = 2

        if (((first_br[1] <= num_r1 and first_br[2] <= num_r1) or (first_br[1] > num_r1 and first_br[2] > num_r1))
            and ra_tot_dict[tet] in first_br[1:3]):

            rel_break = first_br

        elif (((second_br[1] <= num_r1 and second_br[2] <= num_r1) or (second_br[1] > num_r1 and second_br[2] > num_r1))
            and ra_tot_dict[tet] in second_br[1:3]):

            rel_break = second_br

        cl_break = [rel_break[0]]
        for x in rel_break[1:3]:
            if x is not ra_tot_dict[tet]:
                cl_break.append(x)

    # Actual alignment of the molecules after all the preparational work with the rules
    if len(cl_break) != 2:
        break_vec = -rvec_dict[tet][0] # just take one of the tet. reaction vectors as break vector (randomly the first)
    else:
        break_vec = norm_vec(pmg_dict[tet][ra_dict[tet] - 1].coords - \
                             pmg_dict[tet][cl_break[1] - cl_break[0] * num_r1 - 1].coords)

    for vec in rvec_dict[tet]:
        vec_angle = np.absolute(angle_func(vec, break_vec))
        if -5 <= vec_angle <= 5:
            break
        elif 175 <= vec_angle <= 185:
            direction = -1
            break
    rvec_tet = -vec

    return_mol = None

    for rvec_oth in rvec_dict[oth]:

        moltemp_tet = pmg_dict[tet].copy()
        moltemp_oth = pmg_dict[oth].copy()

        R = rot_mat(rvec_oth, rvec_tet)
        z = pmg.core.operations.SymmOp.from_rotation_and_translation(rotation_matrix=R)
        moltemp_oth.apply_operation(z)

        break_vec_rot = direction * norm_vec(moltemp_tet[ra_dict[tet] - 1].coords -
                                             moltemp_tet[cl_break[1] - cl_break[0] * num_r1 - 1].coords)

        t1 = pmg.core.operations.SymmOp.from_rotation_and_translation(
            translation_vec=user_dist * (sum_rads) * break_vec_rot + np.array([0, 0, 1E-3]))
        moltemp_oth.apply_operation(t1)

        if not invert:
            for i in moltemp_oth:
                moltemp_tet.append(species=i.species_string, coords=i.coords, validate_proximity=False)
            return_mol = moltemp_tet
        else:
            for i in moltemp_tet:
                moltemp_oth.append(species=i.species_string, coords=i.coords, validate_proximity=False)
            return_mol = moltemp_oth

    return return_mol
