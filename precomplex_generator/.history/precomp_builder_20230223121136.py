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
import os
import sys
import shutil
import networkx as nx
from openbabel import openbabel as ob
import itertools
from precompiled_fortran_code.get_all_torsions_mod import get_tors
import numpy as np
from tools import angle_func, RingSearch, get_key, norm_vec
from prep_and_align import rvec_align, addSolv
import precompiled_fortran_code.final

from copy import deepcopy
from tools import file_read_dict, file_read_list, file_write, input_gen, reg_precomp, gen_sadd, check_for_rings
from element_dict import element_dict


from platform import python_version

if python_version() < "3.0.0":
    class FileNotFoundError(OSError):
        pass


class preCompBuild():

    """
    Class which contains all information needed for precomplex generation
    """

    def __init__(self, educt_names=None, educts=None, structures=None, topos=None, prec_settings=None, solvent_pmg=None):

        self.structures = structures #dictionary containing pymatgen and openbabel molecule object (and more info)
        self.educt_names = educt_names
        self.educts = educts
        self.topos = topos
        self.prec_settings = prec_settings
        self.solvent_pmg = solvent_pmg


    def mod_conn(self, rule, id, fileIn, fileOut):

        ### Write Connectivity of rule into file

        if len(rule) == 4 or len(rule) == 8:
            ra1 = rule[2] - 1
            ra2 = rule[3] + id - 1
        elif len(rule) == 2:
            ra1 = rule[0] - 1
            ra2 = rule[1] + id - 1

        with open(fileIn, "r") as myfile:
            lines = myfile.readlines()
            temp = []
            for line in lines:
                temp.append(line.split())
            if (str(ra1+1) not in temp[ra2][3:]) and (str(ra2+1) not in temp[ra1][3:]):
                a = int(temp[(ra1)][2]) + 1
                temp[ra1][2] = str(a)
                temp[ra1].append(str(ra2 + 1))
                int_list = list(map(int, temp[ra1][3:]))
                int_list.sort()
                str_list = list(map(str, int_list))
                for i in range(len(temp[ra1][3:])):
                    temp[ra1][i+3] = str_list[i]

                # Modify second atom

                b = int(temp[ra2][2])+1
                temp[ra2][2] = str(b)
                temp[ra2].append(str(ra1 + 1))
                int_list = list(map(int, temp[(ra2)][3:]))
                int_list.sort()
                str_list = list(map(str, int_list))
                for i in range(len(temp[ra2][3:])):
                    temp[ra2][i+3] = str_list[i]
            else:
                print('ERROR: Atoms already connected! Will leave now.')
                pass

            with open(fileOut, "w+") as conn2:
                for lists in temp:
                    conn2.write(' ')
                    for items in lists:
                        conn2.write(items)
                        conn2.write(' ')
                    conn2.write('\n')
            conn2.close()

        return

    def get_disc_tors(self):

        """
        Generate list of torsional angles which will enter conformational search.
        """
        educt_names = self.educt_names
        Structures = self.structures
        disc_tors = []

        for num, educt in enumerate(educt_names):
            if os.path.isfile('../../' + educt + '_tors.dat'):
                temp = file_read_list('../../' + educt + '_tors.dat')
            elif os.path.isfile(educt + '_tors.dat'):
                temp = file_read_list(educt + '_tors.dat')
            else:
                sys.exit("ERROR: _tors.dat cannot be found!")
            for i in range(len(temp)):
                temp[i] = list(map(int, temp[i]))
                if num == 1:
                    temp[i][0] += Structures[0]["Number_of_Atoms"]
                    temp[i][1] += Structures[0]["Number_of_Atoms"]
            disc_tors.extend(temp)

        for i in range(len(disc_tors)):
            disc_tors[i] = list(map(str, disc_tors[i]))

        return disc_tors

    def getNeighborBonds(self, atoms=None, torsList=None, connDict=None, num_spheres=None):

        """
        Go through dictionary of (overall) connectivities (connDict) and only take those which are needed for the
        specific precomplex, because they include bonds to neighbours of the reactive atoms
        (depending on the "number of spheres"=num_spheres (a sphere here is defined as a number of neighbors "distance"
        to be considered, etc.)
        """
        ra_temp = []
        dict_out = {}

        for i in range(int(len(atoms)/2)):
            ra_temp.append(sorted([atoms[i], atoms[i+1]]))

        already = []

        for sphere, i in enumerate(range(num_spheres)):
            registered = []
            new_atoms = []
            for ra in atoms:
                for at in connDict[str(ra)]:
                    if at not in new_atoms and at not in atoms:
                        new_atoms.append(at)
                for single_tors in torsList:
                    if (int(single_tors[0]) == ra or int(single_tors[1]) == ra) and single_tors not in registered \
                            and sorted(list(map(int, single_tors[0:2]))) not in ra_temp and single_tors not in already:
                        registered.append(single_tors)
                        already.append(single_tors)

            atoms = deepcopy(new_atoms)

            dict_out[str(sphere+1)] = registered

        return dict_out

    def modify_conn(self, rule=None, id=None, fileIn=None, fileOut=None):

        """
        Read in connectivity file (fileIn), modify onnectivity and write it again into a new file (fileOut).
        ra1: first reactive atom (of an "add" reaction)
        ra2: second reactive atom (of an "add" reaction)
        """

        if len(rule) == 4 or len(rule) == 8:
            ra1 = rule[2] - 1
            ra2 = rule[3] + id - 1
        elif len(rule) == 2:
            ra1 = rule[0] - 1
            ra2 = rule[1] + id - 1

        with open(fileIn, "r") as myfile:
            lines = myfile.readlines()
            temp = []
            for line in lines:
                temp.append(line.split())
            if (str(ra1+1) not in temp[ra2][3:]) and (str(ra2+1) not in temp[ra1][3:]):
                a = int(temp[(ra1)][2]) + 1
                temp[ra1][2] = str(a)
                temp[ra1].append(str(ra2 + 1))
                int_list = list(map(int, temp[ra1][3:]))
                int_list.sort()
                str_list = list(map(str, int_list))
                for i in range(len(temp[ra1][3:])):
                    temp[ra1][i+3] = str_list[i]

                # Modify second atom
                b = int(temp[ra2][2])+1
                temp[ra2][2] = str(b)
                temp[ra2].append(str(ra1 + 1))
                int_list = list(map(int, temp[(ra2)][3:]))
                int_list.sort()
                str_list = list(map(str, int_list))
                for i in range(len(temp[ra2][3:])):
                    temp[ra2][i+3] = str_list[i]
            else:
                print('ERROR: Atoms already connected! Will leave now.')
                pass

            with open(fileOut, "w+") as conn2:
                for lists in temp:
                    conn2.write(' ')
                    for items in lists:
                        conn2.write(items)
                        conn2.write(' ')
                    conn2.write('\n')
            conn2.close()

        return

    def construct_graph_from_conn(self, conn, ras):

        """
        Construct a graph based on the connectivity information and find shortest-path between two reactive centers.
        """

        edges = []
        for i in conn:
            for j in conn[i]:
                edges.append(sorted([int(i), j]))
                edges.sort()
        edge_set = list(edges for edges, _ in itertools.groupby(edges))
        g = nx.Graph()
        g.add_edges_from(edge_set)

        if len(ras) == 2:
            ra1, ra2 = ras
            path = nx.shortest_path(g, ra1, ra2)
        elif len(ras) == 4:
            ra1, ra2, ra3, ra4 = ras
            path = nx.shortest_path(g, ra1, ra2) + nx.shortest_path(g, ra3, ra4)
        else:
            print("ERROR: NOT implemented! Number of reactive atoms if not equals 2 or 4!")

        path = list(map(int, path))

        return path

    def check_dblelin(self, ras=None, precomp=None, case=None, neighborBonds={}, path=[]):

        """
        Check for planar bonds (=double bonds)
        """
        if case == 'neigh':
            clean_nb = {}

            for dict_key in neighborBonds.keys():
                sphere = neighborBonds[dict_key]
                sphere_temp = []

                for run_id, bond in enumerate(sphere):
                    atom_ids = [int(bond[0]), int(bond[1])]
                    dble_lin = False

                    if ras[0] in atom_ids or ras[1] in atom_ids:
                        set_list = list(set(ras + atom_ids))

                        if len(set_list) == 3:
                            vec1 = precomp[set_list[0]-1].coords - precomp[set_list[1]-1].coords
                            vec2 = precomp[set_list[0]-1].coords - precomp[set_list[2]-1].coords
                            angle = angle_func(vec1, vec2)
                            if abs(angle) < 5 or 175 < abs(angle) < 185:
                                dble_lin = True

                    if not dble_lin:
                        sphere_temp.append(bond)

                clean_nb[dict_key] = sphere_temp

            return clean_nb

        elif case == 'path':
            clean_path = []
            for run_id, bond in enumerate(path):
                atom_ids = [int(bond[0]), int(bond[1])]
                dble_lin = False
                if ras[0] in atom_ids or ras[1] in atom_ids:
                    set_list = list(set(ras + atom_ids))
                    if len(set_list) == 3:
                        vec1 = precomp[set_list[0] - 1].coords - precomp[set_list[1] - 1].coords
                        vec2 = precomp[set_list[0] - 1].coords - precomp[set_list[2] - 1].coords
                        angle = angle_func(vec1, vec2)
                        if abs(angle) < 5 or 175 < abs(angle) < 185:
                            dble_lin = True

                if not dble_lin:
                    clean_path.append(bond)

            return clean_path

    def bond_to_tors(self, bonds=None, tors=[], base_rot=None, planar_fix=True):

        """
        Translate bond/connectivity to torsional angles used for conformational search.
        [planar_fix: if True, DB are not rotated.]

        """

        spheres = list(bonds.keys())
        base_copy = base_rot
        base_copy2 = 45

        if len(bonds) > 0:
            for sphere in spheres:
                factor = int(sphere) # * 2
                for bond in bonds[sphere]:

                    if factor * base_copy < 120:
                        if int(bond[2]) == 120:
                            if [bond[0], bond[1], factor * base_copy] not in tors:
                                tors.append([bond[0], bond[1], factor * base_copy])

                    elif factor * base_copy >= 120:
                        if int(bond[2]) == 120:
                            if [bond[0], bond[1], 120] not in tors:
                                tors.append([bond[0], bond[1], 120])

                    if factor * base_copy2 < 180:
                        if int(bond[2]) == 180:
                            #if planar_fix:
                            #    tors.append([bond[0], bond[1], 180])
                            if not planar_fix:
                                if [bond[0], bond[1], factor * base_copy2] not in tors:
                                    tors.append([bond[0], bond[1], factor * base_copy2])

                    elif factor * base_copy2 >= 180:
                        if int(bond[2]) == 180:
                            if not planar_fix:
                                if [bond[0], bond[1], 180] not in tors:
                                    tors.append([bond[0], bond[1], 180])

        return tors

    def gen_path(self, connDict=None, torsList=None, ras=None):

        path = []
        pre_path = self.construct_graph_from_conn(conn=connDict, ras=ras)

        for i in pre_path:
            for j in torsList:
                if int(j[0]) == i and int(j[1]) in pre_path:
                    path.append(j)

        sort_path = path.copy()

        for i in range(len(sort_path)):
            if i == 0:
                sort_path[i] = path[i]
            elif i % 2 == 0:
                sort_path[i] = path[int(i / 2)]
            elif i % 2 == 1:
                sort_path[i] = path[len(path) - 1 - int(i / 2)]

        return sort_path

    def path_to_tors(self, sort_path, base_rot, tors=[], planar_fix=True):

        """
        Collect torsional dihedral (of a specific value) for the required atom pairs on the sort_path
        """

        base_copy = base_rot
        base_copy2 = 45

        for j, path in enumerate(sort_path):
            if j in [0, 1]:
                if int(path[2]) == 120:
                    tors.append([path[0], path[1], base_copy])
                elif int(path[2]) == 180:
                    #if planar_fix:
                    #    tors.append([path[0], path[1], 180])
                    if not planar_fix:
                        tors.append([path[0], path[1], base_copy2])

            elif j in [2, 3]:
                if int(path[2]) == 120:
                    tors.append([path[0], path[1], 2 * base_copy])
                elif int(path[2]) == 180:
                    #if planar_fix:
                    #    tors.append([path[0], path[1], 180])
                    if not planar_fix:
                        tors.append([path[0], path[1], 2 * base_copy2])

            else:
                if int(path[2]) == 120:
                    tors.append([path[0], path[1], 120])
                elif int(path[2]) == 180:
                    tors.append([path[0], path[1], 180])

        return tors

    def gen_breaks(self, educt_names, topos, breakRule, ras):

        breaks = []
        br_keys = []
        ra_keys = []
        outbreak1 = [0, 0]
        outbreak2 = []

        if breakRule:
            for i in range(int(len(breakRule) / 3)):
                breaks.append([breakRule[3*i], breakRule[3*i+1], breakRule[3*i+2]])

        if len(ras) == 2:
            ra1, ra2 = ras
        elif len(ras) == 4:
            ra1, ra2 = ras[2:4]

        if len(breaks) == 0:
            pass

        else:
            for educt_name in educt_names:
                topo_keys = topos[educt_name].keys()

                for topo_key in topo_keys:

                    if len(breaks) == 1:
                        if (int(topo_key[:3]) == ra1 or int(topo_key[:3]) == ra2) and topo_key[4:7] == 'tet':
                            if ra1 or ra2 in breaks[0]:
                                outbreak1 = [breaks[0][1], breaks[0][2]]

                    if len(breaks) == 2:

                        if int(topo_key[0:3]) == ras[0] or int(topo_key[0:3]) == ras[1]:
                            ra_keys.append(topo_key)

                        elif int(topo_key[0:3]) == breaks[0][1] or int(topo_key[0:3]) == breaks[0][2] or \
                             int(topo_key[0:3]) == breaks[1][1] or int(topo_key[0:3]) == breaks[1][2]:
                             br_keys.append(topo_key)

                        br_keys_set = set(br_keys)
                        ra_keys_set = set(ra_keys)
                        double_at = list(ra_keys_set.intersection(br_keys_set))

                        if len(double_at) > 0 and double_at[0][4:7] == 'tet':
                            if int(double_at[0][:3]) in breaks[0]:
                                outbreak1 = [breaks[0][1], breaks[0][2]]
                            elif int(double_at[0][:3]) in breaks[1]:
                                outbreak1 = [breaks[1][1], breaks[1][2]]

                        else:
                            if int(topo_key[:3]) == ra1 and topo_key[4:7] == 'tet':
                                for sbreak in breaks:
                                    if ra1 in sbreak:
                                        outbreak1 = [sbreak[1], sbreak[2]]

                        if int(topo_key[:3]) == ra2 and topo_key[4:7] == 'tet':
                            for sbreak in breaks:
                                if ra2 in sbreak:
                                    outbreak2 = [sbreak[1], sbreak[2]]

        return [outbreak1, outbreak2]

    def modify_torsList(self, torsList , disc_tors, fadd):


        topos = self.topos
        first_topo = topos[self.educt_names[0]]
        secon_topo = topos[self.educt_names[1]]
        full_topo = []
        for topo_key in first_topo:
            full_topo.append(topo_key)
        for topo_key in secon_topo:
            full_topo.append(topo_key)

        dis_con_ids = []
        dis_con_angles = []
        for i in disc_tors:
            dis_con_ids.append([i[0], i[1]])
            dis_con_angles.extend(i[2])

        torsList_ids = []
        for i in torsList:
            torsList_ids.append([i[0], i[1]])
            dis_con_angles.extend(i[2])

        final_torsList = []
        for tors in disc_tors:
            if ((full_topo[int(tors[0])-1][4:7] == 'tpl') and (full_topo[int(tors[1])-1][4:7] == 'tpl')):
                final_torsList.append([tors[0], tors[1], '180'])
            else:
                final_torsList.append(tors)

        for tors in torsList:
            if sorted([tors[0], tors[1]]) == sorted(list(map(str, fadd))):
                final_torsList.append(tors)


        diff_tors = abs(len(final_torsList) - len(torsList))
        diff_ids = abs(len(final_torsList) - len(torsList_ids))

        if diff_tors == 0 and diff_ids == 0:
            pass

        elif diff_tors > 0 and diff_ids > 0:
            miss_tors = []
            for i in torsList:
                if i not in final_torsList:
                    miss_tors.append(i)

            for miss_tor in miss_tors:

                if (((full_topo[int(miss_tor[0])-1][4:7]) == 'tpl') and ((full_topo[int(miss_tor[1])-1][4:7]) == 'tpl')):
                    final_torsList.append([miss_tor[0], miss_tor[1], '180'])
                    # print('New pre-complex bond on planar system. Will still rotate around 180°!')

                else:
                    final_torsList.append(miss_tor)
        else:
            print('Error in modify_torsList: Case not implemented.')

        return final_torsList


    def intramolReaction(self, topos=None, addRule=None, breakRule=None, xyzFile=None, preCompList=None):

        """
        Prepare precomplex generation for intramolecular reactions.
        """

        intramol_small_ring = False
        base_rot = self.prec_settings['base_rot']
        num_spheres = self.prec_settings['num_spher']

        fortran_weights = np.asarray([
            self.prec_settings['max_min_weight'],
            self.prec_settings['avg_dist_weight'],
            self.prec_settings['add_dist_weight'],
            self.prec_settings['angle_diff_weight'],
            self.prec_settings['dummy_weight'],
            self.prec_settings['cylin_weight'],
            self.prec_settings['numrot_weight'],],
            dtype=np.float32)

        sadd_factor = self.prec_settings['sadd_weight']
        planar_fix = self.prec_settings['planar_fix']

        max_mem = 10
        iter = 0
        educt_names = self.educt_names

        connDict = file_read_dict(educt_names[0] + '_conn.dat')
        get_tors(educt_names[0] + '_conn.dat')
        os.rename('tors.dat', educt_names[0] + '_tors.dat')
        torsList = file_read_list(educt_names[0] + '_tors.dat')

        xyzAppend = []

        # case of 1-ADD
        if len(addRule) == 4:

            ras = list(filter((0).__ne__, addRule))
            ringsize, rings = RingSearch(obmols=self.educts, atIDs=[[ras[0], ras[1]]])
            if ringsize <= 6: # Note: 6 might not always be the best criterion here!
                intramol_small_ring = True

            topo1_keys = topos[educt_names[0]].keys()
            try:
                key1 = get_key(inp_dict=topo1_keys, ra=ras[0])
                key2 = get_key(inp_dict=topo1_keys, ra=ras[1])

                if ringsize <= 6:
                    rvecs1 = topos[educt_names[0]][key1]["small_rings"]
                    rvecs2 = topos[educt_names[0]][key2]["small_rings"]
                else:
                    rvecs1 = topos[educt_names[0]][key1]["normal"]
                    rvecs2 = topos[educt_names[0]][key2]["normal"]
            except:
                print("Problem with Precomplex generation...!")
                return preCompList

            breakvec = None
            rvec1_needed = False
            rvec2_needed = False

            # check for SN2-like add and break situation; i.e. one of the add atoms must occur in one of the break atom lists!
            if breakRule != None and breakRule !=[]:
                if len(breakRule) == 3:
                    # check if and which "add atom" occurs in "breaks"
                    if addRule[2] == breakRule[1]:
                        break_a1 = breakRule[1]
                        break_a2 = breakRule[2]
                        rvec1_needed = True
                    elif addRule[2] == breakRule[2]:
                        break_a1 = breakRule[2]
                        break_a2 = breakRule[1]
                        rvec1_needed = True
                    elif addRule[3] == breakRule[1]:
                        break_a1 = breakRule[1]
                        break_a2 = breakRule[2]
                        rvec2_needed = True
                    elif addRule[3] == breakRule[2]:
                        break_a1 = breakRule[2]
                        break_a2 = breakRule[1]
                        rvec2_needed = True
                elif len(breakRule) == 6:
                    if addRule[2] == breakRule[1]:
                        break_a1 = breakRule[1]
                        break_a2 = breakRule[2]
                        rvec1_needed = True
                    elif addRule[2] == breakRule[2]:
                        break_a1 = breakRule[2]
                        break_a2 = breakRule[1]
                        rvec1_needed = True
                    elif addRule[3] == breakRule[1]:
                        break_a1 = breakRule[1]
                        break_a2 = breakRule[2]
                        rvec2_needed = True
                    elif addRule[3] == breakRule[2]:
                        break_a1 = breakRule[2]
                        break_a2 = breakRule[1]
                        rvec2_needed = True
                    elif addRule[2] == breakRule[4]:
                        break_a1 = breakRule[4]
                        break_a2 = breakRule[5]
                        rvec1_needed = True
                    elif addRule[2] == breakRule[5]:
                        break_a1 = breakRule[5]
                        break_a2 = breakRule[4]
                        rvec1_needed = True
                    elif addRule[3] == breakRule[4]:
                        break_a1 = breakRule[4]
                        break_a2 = breakRule[5]
                        rvec2_needed = True
                    elif addRule[3] == breakRule[5]:
                        break_a1 = breakRule[4]
                        break_a2 = breakRule[4]
                        rvec2_needed = True

                xyz_list = file_read_list(filename=self.educt_names[0] + '.xyz')
                atomCoord1 = np.asarray(list(map(float, xyz_list[int(break_a1) - 1][1:])))
                atomCoord2 = np.asarray(list(map(float, xyz_list[int(break_a2) - 1][1:])))

                breakvec = [norm_vec(atomCoord2 - atomCoord1), norm_vec(atomCoord2 - atomCoord1),
                            norm_vec(atomCoord2 - atomCoord1), norm_vec(atomCoord2 - atomCoord1)]


            for i in range(len(rvecs1)):
                for j in range(len(rvecs2)):


                    if breakvec != None:
                        if rvec1_needed:
                            crossp_norm = np.linalg.norm(np.cross(rvecs1[i], breakvec))
                        elif rvec2_needed:
                            crossp_norm = np.linalg.norm(np.cross(rvecs2[i], breakvec))

                        # for SN2-like add/break situation, only use the vector of the 4 "tetrahedral" ones, which points
                        # away from the bond to be broken!
                        if crossp_norm != 0:
                            continue

                    geomFolder = 'geom_' + str(iter)
                    iter += 1

                    if not os.path.exists(geomFolder):
                        os.makedirs(geomFolder)

                    os.chdir(geomFolder)

                    shutil.copyfile('../' + educt_names[0] + '_conn.dat', './' + educt_names[0] + '_conn.dat')
                    shutil.copyfile('../' + educt_names[0] + '_tors.dat', './' + educt_names[0] + '_tors.dat')
                    shutil.copyfile('../' + educt_names[0] + '.xyz', './' + educt_names[0] + '.xyz')

                    path = self.gen_path(connDict=connDict, torsList=torsList, ras=ras)

                    input_gen(new_data=[ras[0], ras[1], gen_sadd(path=path, sadd_factor=sadd_factor)],
                              data=xyzAppend, case='add')

                    out_breaks = self.gen_breaks(educt_names=educt_names, topos=topos, breakRule=breakRule, ras=ras)
                    input_gen(new_data=out_breaks, data=xyzAppend, case='breaks')

                    tors = self.path_to_tors(sort_path=path, base_rot=base_rot, planar_fix=planar_fix)

                    neighborBonds = self.getNeighborBonds(atoms=ras, torsList=torsList, num_spheres=num_spheres,
                                                          connDict=connDict)
                    self.bond_to_tors(bonds=neighborBonds, tors=tors, base_rot=base_rot, planar_fix=planar_fix)

                    input_gen(new_data=tors, data=xyzAppend, case='tors')

                    file_write(final_file=xyzFile, final_data=xyzAppend, case='intramol')

                    final.file_input(xyzFile, educt_names[0] + '_conn.dat', 'intramol', fortran_weights, max_mem,
                                     len(neighborBonds), len(path), i, j, intramol_small_ring)

                    reg_precomp(preCompList=preCompList, input_name=educt_names[0])

                    os.chdir('..')

        # case of 2-ADDs
        elif len(addRule) == 8:

            ras_tmp = list(filter((0).__ne__, addRule))
            topo1_keys = topos[educt_names[0]].keys()

            # currently: randomly choose one add as the first one
            # better: use both orderings and take best overall precomplex structure!

            ras = [ras_tmp[0], ras_tmp[1]]

            ringsize, rings = RingSearch(obmols=self.educts, atIDs=[[ras_tmp[0], ras_tmp[1], ras_tmp[2], ras_tmp[3]]])

            if len(rings) > 1: # more than 1 ring formed!

                outside = False

                if len(rings[0]) == len(rings[1]): # it cannot be that this is one ring within another...
                    outside = True

                if len(rings[0]) < len(rings[1]):
                    small_ring = rings[0]
                    for it in rings[0]:
                        if it not in rings[1]:
                            outside = True

                elif len(rings[1]) < len(rings[0]):
                    small_ring = rings[1]
                    for it in rings[1]:
                        if it not in rings[0]:
                            outside = True


                if not outside: # the two rings "overlap" regarding the atomIDs in between the reactive atoms
                    ras = [small_ring[0], small_ring[-1]]

            try:
                key1 = get_key(inp_dict=topo1_keys, ra=ras[0])
                key2 = get_key(inp_dict=topo1_keys, ra=ras[1])
                rvecs1 = topos[educt_names[0]][key1]
                rvecs2 = topos[educt_names[0]][key2]
            except:
                print("Problem with Precomplex generation...!")
                return preCompList

            for i in range(len(rvecs1)):
                for j in range(len(rvecs2)):

                    geomFolder = 'geom_' + str(iter)
                    iter += 1

                    if not os.path.exists(geomFolder):
                        os.makedirs(geomFolder)

                    os.chdir(geomFolder)

                    shutil.copyfile('../' + educt_names[0] + '_conn.dat', './' + educt_names[0] + '_conn.dat')
                    shutil.copyfile('../' + educt_names[0] + '_tors.dat', './' + educt_names[0] + '_tors.dat')
                    shutil.copyfile('../' + educt_names[0] + '.xyz', './' + educt_names[0] + '.xyz')

                    path = self.gen_path(connDict=connDict, torsList=torsList, ras=ras)

                    input_gen(new_data=[ras[0], ras[1], gen_sadd(path=path, sadd_factor=sadd_factor)],
                              data=xyzAppend, case='add')
                    out_breaks = self.gen_breaks(educt_names=educt_names, topos=topos, breakRule=breakRule, ras=ras)
                    input_gen(new_data=out_breaks, data=xyzAppend, case='breaks')
                    self.path_to_tors(sort_path=path, tors=tors, base_rot=base_rot, planar_fix=planar_fix)
                    neighborBonds = self.getNeighborBonds(atoms=ras, torsList=torsList, num_spheres=num_spheres,
                                                              connDict=connDict)
                    self.bond_to_tors(bonds=neighborBonds, tors=tors, base_rot=base_rot, planar_fix=planar_fix)

                    input_gen(new_data=tors, data=xyzAppend, case='tors')

                    file_write(final_file=xyzFile, final_data=xyzAppend, case='intramol')

                    final.file_input(xyzFile, educt_names[0] + '_conn.dat', 'intramol', fortran_weights, max_mem,
                                     len(neighborBonds), len(path), i, j, intramol_small_ring)

                    reg_precomp(preCompList=preCompList, input_name=educt_names[0])

                    os.chdir('..')

        elif len(addRule) == 0:
            reg_precomp(preCompList=preCompList, input_name=educt_names[0])
        else:
            reg_precomp(preCompList=preCompList, input_name=educt_names[0])

        return preCompList


    def bimolecularReaction(self, educt_names=None, structures=None, topos=None, addRule=None, breakRule=None,
                            xyzFile=None, preCompList=None):


        """
        Prepare precomplex generation for bimolecular reactions.
        """

        intramol_small_ring = False

        if len(addRule) < 4:
            return None

        ras = list(filter((0).__ne__, addRule))

        # find out if second add has "prochiral" centers (i.e. len(rvecs) > 1)...
        if len(addRule) > 4:
            topo1_keys = topos[educt_names[addRule[4]]].keys()
            topo2_keys = topos[educt_names[addRule[5]]].keys()
            try:
                key1 = get_key(inp_dict=topo1_keys, ra=addRule[6])
                key2 = get_key(inp_dict=topo2_keys, ra=addRule[7])
                rvecs1 = topos[educt_names[addRule[4]]][key1]["normal"] # educt name id of second add...
                rvecs2 = topos[educt_names[addRule[5]]][key2]["normal"] # educt name id of second add...
            except:
                print("Problem with Precomplex generation...!")
                return preCompList
        else:
            rvecs1 = []
            rvecs2 = []

        base_rot = self.prec_settings['base_rot']
        num_spheres = self.prec_settings['num_spher']
        dist_weight = self.prec_settings['dist_weight']

        nTotal = structures[0]["Number_of_Atoms"] + structures[1]["Number_of_Atoms"]
        self.modify_conn(rule=addRule, id=structures[0]["Number_of_Atoms"], fileIn="dis_conn.dat",
                     fileOut="precomp_conn_fa.dat")

        ### Comparison of torsional angles for reactants and precomplex

        get_tors('precomp_conn_fa.dat')
        os.rename('tors.dat', 'conn_tors.dat')

        connDict = file_read_dict("precomp_conn_fa.dat")
        torsList = file_read_list('conn_tors.dat')

        disc_tors = self.get_disc_tors()

        torsList = self.modify_torsList(torsList=torsList, disc_tors=disc_tors,
                                        fadd=[addRule[2] + addRule[0] * structures[0]["Number_of_Atoms"],
                                              addRule[3] + addRule[1] * structures[0]["Number_of_Atoms"]])

        intelligentDistance = dist_weight * (2 + 0.5*np.log10(nTotal / 10))   # Note: This threshold might need to be adapted!

        geomBiMol, used_rvecs, topo_keys = rvec_align(educt_names=educt_names, educts=structures, addRule=addRule[0:4], breakRules=breakRule,
                                                      topos=topos, user_dist=intelligentDistance)

        fortran_weights = np.asarray([
            self.prec_settings['max_min_weight'],
            self.prec_settings['avg_dist_weight'],
            self.prec_settings['add_dist_weight'],
            self.prec_settings['angle_diff_weight'],
            self.prec_settings['dummy_weight'],
            self.prec_settings['cylin_weight'],
            self.prec_settings['numrot_weight'],],
            dtype=np.float32)

        sadd_factor = self.prec_settings['sadd_weight']
        planar_fix = self.prec_settings['planar_fix']

        max_mem = 20
        iter = 0

        if len(addRule) == 4:

            for i in range(len(geomBiMol)):

                xyzAppend = []
                geomFolder = 'geom_' + str(i)

                if not os.path.exists(geomFolder):
                    os.makedirs(geomFolder)
                os.chdir(geomFolder)

                if geomBiMol[i] == None:
                    continue
                try:
                    geomBiMol[i].to(fmt='xyz', filename='precomp_pre_conf.xyz')
                except AttributeError:
                    pass
                shutil.copyfile('../precomp_conn_fa.dat', './precomp_conn_fa.dat')

                ras = [addRule[2], addRule[3] + structures[0]["Number_of_Atoms"]]
                input_gen(new_data=[ras[0], ras[1], base_rot], data=xyzAppend, case='fadd')

                out_breaks = self.gen_breaks(educt_names=educt_names, topos=topos, breakRule=breakRule, ras=ras)
                input_gen(new_data=out_breaks, data=xyzAppend, case='breaks')

                neighborBonds_pre = self.getNeighborBonds(atoms=ras, torsList=torsList, num_spheres=num_spheres,
                                                          connDict=connDict)
                neighborBonds = self.check_dblelin(ras=ras, neighborBonds=neighborBonds_pre,
                                                   precomp=geomBiMol[i], case='neigh')

                tors = self.bond_to_tors(bonds=neighborBonds, base_rot=base_rot, planar_fix=planar_fix)

                input_gen(new_data=tors, data=xyzAppend, case='tors')
                file_write(final_file=xyzFile, final_data=xyzAppend, case='oneadd')
                final.file_input('precomp_pre_conf.xyz', 'precomp_conn_fa.dat', 'oneadd', fortran_weights, max_mem,
                                 len(neighborBonds), 0, 0, 0, intramol_small_ring)

                reg_precomp(preCompList=preCompList)
                os.chdir('..')

        elif len(addRule) == 8:

            rule = {addRule[0]:[], addRule[1]:[]}
            rule[addRule[0]].append(addRule[2])
            rule[addRule[1]].append(addRule[3])
            rule[addRule[4]].append(addRule[6])
            rule[addRule[5]].append(addRule[7])

            # check if reactive atoms of first/second educt molecule are in a ring!
            ring_check1, neighbors1 = check_for_rings(educt=self.educts[addRule[0]]["obm"], at_num_1=rule[0][0], at_num_2=rule[0][1])
            ring_check2, neighbors2 = check_for_rings(educt=self.educts[addRule[1]]["obm"], at_num_1=rule[1][0], at_num_2=rule[1][1])

            iterator = 0

            for i in range(len(geomBiMol)):
                for a in range(len(rvecs1)):
                    for b in range(len(rvecs2)):

                        iterator += 1
                        geomFolder = 'geom_' + str(iter)
                        iter += 1

                        xyzAppend = []

                        #################
                        # it must be checked if first add was also a "tpl" center (found in topo_keys)

                        if "tpl" in key1:
                            if "tpl" in topo_keys[0]:
                                for key in used_rvecs[i][0].keys():
                                    len_vec1 = np.linalg.norm(rvecs1[a])
                                    crossp_norm1 = np.linalg.norm(np.cross(rvecs1[a], used_rvecs[i][0][key]))
                                    sum1 = np.linalg.norm(np.add(rvecs1[a], used_rvecs[i][0][key]))

                        if "tpl" in key2:
                            if "tpl" in topo_keys[1]:
                                for key in used_rvecs[i][1].keys():
                                    len_vec2 = np.linalg.norm(rvecs2[b])
                                    crossp_norm2 = np.linalg.norm(np.cross(rvecs2[b], used_rvecs[i][1][key]))
                                    sum2 = np.linalg.norm(np.add(rvecs2[b],used_rvecs[i][1][key]))

                        # Do not generate precomplexes with reaction vectors in antiparallel direction, if the tpl reaction centers are
                        # A) neighbor atoms (e.g. C=C double bond) or B) are part of a ring structure
                        # regarding B: one would need to check the ring size, to be sure that an antiparallel reaction vector alignment can
                        # safely be neglected...

                        if "tpl" in key1 and (ring_check1 or neighbors1):
                            if crossp_norm1 < 0.4 and sum1 < len_vec1:
                                continue

                        if "tpl" in key2 and (ring_check2 or neighbors2):
                            if crossp_norm2 < 0.4 and sum2 < len_vec2:
                                continue

                        #################
                        if geomBiMol[i] == None:
                            continue

                        #geomFolder = 'geom_' + str(i)
                        if not os.path.exists(geomFolder):
                            os.makedirs(geomFolder)
                        os.chdir(geomFolder)

                        try:
                            geomBiMol[i].to(fmt='xyz', filename='precomp_pre_conf.xyz')
                        except AttributeError:
                            pass

                        shutil.copyfile('../precomp_conn_fa.dat', './precomp_conn_fa.dat')

                        ras = [addRule[2] + addRule[0] * structures[0]["Number_of_Atoms"], addRule[3] + addRule[1] * structures[0]["Number_of_Atoms"],
                               addRule[6] + addRule[4] * structures[0]["Number_of_Atoms"], addRule[7] + addRule[5] * structures[0]["Number_of_Atoms"]]

                        input_gen(new_data=[ras[0], ras[1], base_rot], data=xyzAppend, case='fadd')

                        out_breaks = self.gen_breaks(educt_names=educt_names, topos=topos, breakRule=breakRule, ras=ras)
                        input_gen(new_data=out_breaks, data=xyzAppend, case='breaks')

                        path_pre = self.gen_path(connDict=connDict, torsList=torsList, ras=ras[2:4])
                        path = self.check_dblelin(ras=[ras[0], ras[1]], path=path_pre, precomp=geomBiMol[i], case='path')
                        tors = self.path_to_tors(sort_path=path, base_rot=base_rot)

                        input_gen(new_data=[ras[2], ras[3], gen_sadd(path=path, sadd_factor=sadd_factor)], data=xyzAppend,
                                  case='add')

                        neighborBonds_pre = self.getNeighborBonds(atoms=ras, torsList=torsList, num_spheres=num_spheres,
                                                                  connDict=connDict)

                        neighborBonds = self.check_dblelin(ras=[ras[0], ras[1]], neighborBonds=neighborBonds_pre,
                                                           precomp=geomBiMol[i], case='neigh')

                        self.bond_to_tors(bonds=neighborBonds, base_rot=base_rot, tors=tors, planar_fix=planar_fix)

                        input_gen(new_data=tors, data=xyzAppend, case='tors')

                        file_write(final_file=xyzFile, final_data=xyzAppend, case='twoadds')

                        final.file_input('precomp_pre_conf.xyz', 'precomp_conn_fa.dat', 'twoadds', fortran_weights,
                                          max_mem, len(neighborBonds), len(path), a, b, intramol_small_ring)

                        reg_precomp(preCompList=preCompList)
                        os.chdir('..')

        elif len(addRule) == 0:
            reg_precomp(preCompList=preCompList)

        return preCompList


    def write_out_solvconn(self, educts, shuttle_mol):

        name = "disc_solv.dat"
        out = open(name, "w")
        neighbour_dic = {}
        shift1 = educts[0]["Number_of_Atoms"]
        shift2 = educts[1]["Number_of_Atoms"]
        for i, mol in enumerate(educts):
            indexShift = i * shift1
            for atom in ob.OBMolAtomIter(mol["obm"]):
                id1 = atom.GetIdx() + indexShift
                # symbol1 = ob.OBElementTable().GetSymbol(atom.GetAtomicNum())
                symbol1 = element_dict[atom.GetAtomicNum()]
                tmp = []
                nneigh = 0
                for neigh in ob.OBAtomAtomIter(atom):
                    nneigh += 1
                    id2 = neigh.GetIdx() + indexShift
                    neighbour_string = str(id2)
                    tmp.append(neighbour_string)
                neighbour_dic[str(id1)] = str(id1) + " " + symbol1 + " " + str(nneigh) + " " + " ".join(sorted(tmp))
        for atom in ob.OBMolAtomIter(shuttle_mol["obm"]):
            indexShift = shift1 + shift2
            id1 = atom.GetIdx() + indexShift
            # symbol1 = ob.OBElementTable().GetSymbol(atom.GetAtomicNum())
            symbol1 = element_dict[atom.GetAtomicNum()]
            tmp = []
            nneigh = 0
            for neigh in ob.OBAtomAtomIter(atom):
                nneigh += 1
                id2 = neigh.GetIdx() + indexShift
                neighbour_string = str(id2)
                tmp.append(neighbour_string)
            neighbour_dic[str(id1)] = str(id1) + " " + symbol1 + " " + str(nneigh) + " " + " ".join(sorted(tmp))
        for key in neighbour_dic.keys():
            out.write(neighbour_dic[key] + "\n")
        out.close()

    def solventShuttle(self, educts=None, educt_names=None, structures=None, topos=None, addRule=None, breakRule=None,
                       preCompList=None, shuttle=None, shuttle_mol=None):

        solvent_pmg = self.solvent_pmg

        intramol_small_ring = False

        fortran_weights = np.asarray([
            self.prec_settings['max_min_weight'],
            self.prec_settings['avg_dist_weight'],
            self.prec_settings['add_dist_weight'],
            self.prec_settings['angle_diff_weight'],
            self.prec_settings['dummy_weight'],
            self.prec_settings['cylin_weight'],
            self.prec_settings['numrot_weight'], ],
            dtype=np.float32)

        sadd_factor = self.prec_settings['sadd_weight']

        max_mem = 10

        if len(educts) == 1 or len(educts) == 2:

            base_rot = self.prec_settings['base_rot']
            num_spheres = self.prec_settings['num_spher']
            dist_weight = self.prec_settings['dist_weight']
            tors = []

            if len(educts) == 2:
                nTotal = educts[0]["Number_of_Atoms"] + educts[1]["Number_of_Atoms"]
            else:
                nTotal = educts[0]["Number_of_Atoms"]

            add_1, add_2 = addRule[0:4], addRule[4:8]
            symbs_1, symbs_2 = [], []

            #print("add_1, add_2:", add_1, add_2)

            for i in range(2, 4):
                symbs_1.append(educts[add_1[i-2]]["pmg"][add_1[i]-1].species_string)
                if len(educts) == 2:
                    symbs_2.append(educts[add_2[i-2]]["pmg"][add_2[i]-1].species_string)

            if 'H' in symbs_1 and len(educts) == 2:
                temp = add_1.copy()
                add_1 = add_2.copy()
                add_2 = temp.copy()

            if len(educts) == 1:
                self.mod_conn(rule=add_1, id=0, fileIn=educt_names[0]+"_conn.dat", fileOut="precomp_conn_fa.dat")
            else:
                self.mod_conn(rule=add_1, id=len(structures[0]), fileIn="dis_conn.dat", fileOut="precomp_conn_fa.dat")

            get_tors('precomp_conn_fa.dat')

            os.system('mv tors.dat tors_backup.dat')

            #ToDo: if len(educts) == 1!!!
            if len(educts) == 1:
                print("Warning: No shuttle for unimolecular reactions implemented!")
                print("Submit as bimolecular reaction!")

            elif len(educts) == 2:

                intelligentDistance = dist_weight * (2 + 0.5*np.log10(nTotal / 10))  # Muss eventuell vergrößert werden

                geomBiMol, used_rvecs, topo_keys = rvec_align(educt_names=educt_names, educts=structures, addRule=add_1, breakRules=breakRule,
                                                              topos=topos, user_dist=intelligentDistance)

                disc_geomBiMol, used_rvecs, topo_keys = rvec_align(educt_names=educt_names, educts=structures, addRule=add_1, breakRules=breakRule,
                                                                   topos=topos, user_dist=1000)

                cplx_acc = add_2[2] + add_2[0] * len(structures[0]["pmg"])
                cplx_don = add_2[3] + add_2[1] * len(structures[0]["pmg"])
                solv_acc = int(shuttle[2]) + nTotal
                solv_don = int(shuttle[1]) + nTotal

                for i in range(0, len(geomBiMol)):
                    add_1 = add_1.copy()

                    xyzAppend = []

                    geomFolder = 'geom_' + str(i)

                    if not os.path.exists(geomFolder):
                        os.makedirs(geomFolder)
                    os.chdir(geomFolder)

                    geomBiMol[i].to(fmt='xyz', filename='precomp_pre_solv.xyz')

                    temp_fa_cplx = disc_geomBiMol[i].copy()
                    temp_solv = solvent_pmg.copy()
                    temp_solv.translate_sites(vector=[0, 0, 1E9])

                    for solv_site in temp_solv:
                        temp_fa_cplx.append(species=solv_site.species_string, coords=solv_site.coords, validate_proximity=False)

                    temp_fa_cplx.to(fmt='xyz', filename='disc_solv.xyz')
                    self.write_out_solvconn(educts, shuttle_mol[0])

                    self.mod_conn(rule=add_1, id=len(structures[0]["pmg"]), fileIn="disc_solv.dat",
                                 fileOut="precomp_conn_fa.dat")

                    solv_cplx = addSolv(pmgs=[geomBiMol[i], solvent_pmg], solv_acc=solv_acc,
                                        cplx_don=cplx_don, topos=topos, shuttle_tup=shuttle,
                                        prev_ras=add_1, len_first=len(structures[0]["pmg"]), len_both=nTotal,
                                        intelligentDistance=intelligentDistance)

                    self.mod_conn(rule=[cplx_don, solv_acc-nTotal], id=nTotal, fileIn="precomp_conn_fa.dat",
                                 fileOut="precomp_solv_add.dat")

                    conn_dict = file_read_dict('precomp_solv_add.dat')
                    get_tors('precomp_solv_add.dat')
                    torsList = file_read_list('tors.dat')

                    # Generate new rule with solvent molecule included!
                    # Therefore, modify existing rule!
                    self.solventRule = {"ADD": [[add_1[0], add_1[1], add_1[2] + add_1[0] * len(structures[0]["pmg"]),
                                                add_1[3] + add_1[1] * len(structures[0]["pmg"])],
                                                [add_2[1], 2, int(cplx_don), int(solv_acc)], # add_2[1] = cplx_don_ID!
                                                [add_2[0], 2, int(cplx_acc), int(solv_don)]],
                                        "BREAK": [[2, int(solv_acc), int(solv_don)],
                                                  [breakRule[0], breakRule[1]+breakRule[0]*len(structures[0]["pmg"]),
                                                   breakRule[2]+breakRule[0]*len(structures[0]["pmg"])]]}

                    input_gen(new_data=[add_1[2] + add_1[0] * len(structures[0]["pmg"]), add_1[3] + add_1[1] * len(structures[0]["pmg"]),
                              base_rot], data=xyzAppend, case='fadd')



                    path_pre = self.gen_path(connDict=conn_dict, torsList=torsList, ras=[cplx_don, solv_acc])
                    path = self.check_dblelin(ras=[cplx_don, solv_acc], path=path_pre, precomp=solv_cplx, case='path')
                    self.path_to_tors(sort_path=path, tors=tors, base_rot=base_rot)


                    input_gen(new_data=[cplx_acc, solv_don, gen_sadd(path=path, sadd_factor=sadd_factor)], data=xyzAppend,
                              case='add')

                    neighborBonds_pre = self.getNeighborBonds(atoms=[cplx_don, int(shuttle[2])], torsList=torsList,
                                                              num_spheres=num_spheres, connDict=conn_dict)
                    neighborBonds = self.check_dblelin(ras=[cplx_don, int(shuttle[2])], neighborBonds=neighborBonds_pre,
                                                       precomp=solv_cplx, case='neigh')

                    self.bond_to_tors(bonds=neighborBonds, tors=tors, base_rot=base_rot)

                    input_gen(new_data=tors, data=xyzAppend, case='tors')

                    file_write(final_file='cplx_with_solv.xyz', final_data=xyzAppend, case='shuttle')

                    final.file_input('cplx_with_solv.xyz', 'precomp_solv_add.dat', 'shuttle', fortran_weights,
                                     max_mem, len(neighborBonds), len(path), 0, 0, intramol_small_ring)

                    reg_precomp(preCompList=preCompList)

                    os.chdir('..')

        return preCompList
