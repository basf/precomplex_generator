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

import argparse
from .tools import add_rule_stand, break_rule_stand, getEducts, genDisConn
from .prep_and_align import write_out_disconn, topo_analysis
from .precomp_builder import preCompBuild
import os, sys
import shutil
import json
from pymatgen import Molecule


def prep_topo(educts=None, educt_names=None, own_angle=None, shuttle=None, shuttle_mol=None):

    """
    Tasks:
    1) Analyse the topology (number of neighbors, etc.) for each reactant structure
    2) generate disconnected structures (for bimolecular reactions)
    """

    topoData = {}

    if len(educt_names) == 1:
        topoData[educt_names[0]] = topo_analysis(educt_name=educt_names[0], educt=educts[0], own_angle=own_angle)

    elif len(educt_names) == 2:
        for i, educt_name in enumerate(educt_names):
            topoData[educt_name] = topo_analysis(educt_name=educt_name, educt=educts[i], own_angle=own_angle)

    if shuttle != None:
        if shuttle[0] != "None":
            topoData[shuttle[0]] = topo_analysis(educt_name=shuttle[0], educt=shuttle_mol[0], shuttle=True)

    ### Make class instance of Pre-Complex-Builder class
    genDisConn(structures=educts)
    write_out_disconn(educts)

    return topoData

def generatePreComplex(rule, educt_names=None, educts=None, build_class=None, shuttle=None, shuttle_mol=None):
    """
    Generate Pre-Complexes for uni- and bimolecular reactions.
    """

    workdir = os.getcwd()

    ### Do not remove buildPC folder if exists, create if not
    if not os.path.exists('buildPC'):
        os.makedirs('buildPC')

    # standardize add and break rules
    addRule = add_rule_stand(structures=educts, rule=rule)
    breakRule = break_rule_stand(structures=educts, rule=rule)

    execPath = workdir
    rule_path = os.path.join(workdir, 'buildPC', "PATH_0")

    # if no topology analysis could be done, we cannot proceed further...
    if build_class == None:
        return

    topos = build_class.topos

    # clean up folder, if paths already exist!
    if os.path.exists(rule_path):
        shutil.rmtree(rule_path)
        os.makedirs(rule_path)
    else:
        os.makedirs(rule_path)

    os.chdir(rule_path)

    # unimolecular/intramolecular reaction
    if len(educt_names) == 1 and shuttle_mol == None:
        shutil.copyfile(execPath + '/' + educt_names[0] + '.xyz', rule_path + '/' + educt_names[0] + '.xyz')
        shutil.copyfile(execPath + '/' + educt_names[0] + '_conn.dat',
                        rule_path + '/' + educt_names[0] + '_conn.dat')
        precomplexes = build_class.intramolReaction(topos=topos, addRule=addRule, breakRule=breakRule,
                                                    xyzFile=educt_names[0] + '.xyz', preCompList=[])
    # bimolecular reaction
    elif len(educt_names) == 2 and shuttle_mol == None:
        shutil.copyfile(execPath + '/dis_conn.dat', rule_path + '/dis_conn.dat')
        precomplexes = build_class.bimolecularReaction(educt_names=educt_names, structures=educts, topos=topos,
                                                       addRule=addRule, breakRule=breakRule,
                                                       xyzFile='precomp_pre_conf.xyz', preCompList=[])

    elif len(educt_names) == 1 and shuttle_mol != None:
        shutil.copyfile(execPath + '/' + educt_names[0] + '.xyz', rule_path + '/' + educt_names[0] + '.xyz')
        shutil.copyfile(execPath + '/' + educt_names[0] + '_conn.dat', rule_path + '/' + educt_names[0] + '_conn.dat')

        precomplexes = build_class.solventShuttle(topos=topos, addRule=addRule, breakRule=breakRule,
                                                  educt_names=educt_names, educts=educts,
                                                  preCompList=[], shuttle=shuttle, shuttle_mol=shuttle_mol)

    elif len(educt_names) == 2 and len(addRule) == 8 and shuttle_mol != None:
        shutil.copyfile(execPath + '/dis_conn.dat', rule_path + '/dis_conn.dat')
        precomplexes = build_class.solventShuttle(educts=educts, educt_names=educt_names, structures=educts,
                                                  topos=topos, addRule=addRule, breakRule=breakRule,
                                                  preCompList=[], shuttle=shuttle, shuttle_mol=shuttle_mol)

    os.chdir(workdir)

    return precomplexes


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='QuantumNet: a Demonstrator...',
    )
    parser.add_argument('educts', type=str, nargs='*',
                        help="Reactants: molecule1 and/or molecule2 (name of XYZ files).")
    parser.add_argument('-r', '--rule', dest='r', default=None, type=str,
                        help='Specification of one single rule to be calculated')
    args = parser.parse_args()

    types = ["xyz", "xyz"]
    if not len(args.educts) > 0 or not len(types) > 0:
        raise SystemError("Please provide at least 1 input molecule!")

    educt_list = []

    for ed in args.educts:
       if ed[-4:] == ".xyz":
          educt_list.append(ed[:-4])
       else:
          educt_list.append(ed)

    educt_names = sorted(educt_list)

    mol = [{"name": educt_names[0], "type": types[0]}]
    if len(args.educts) > 1 and len(types) > 1:
        mol.append({"name": educt_names[1], "type": types[1]})

    # get pymatgen and openbabel molecule objects for the provided reactant XYZ structures
    educts = getEducts(mol, verbose=True)

    if len(educts) == 2:
        natoms = educts[0]["Number_of_Atoms"]

    adds = []
    breaks = []
    subst_names = []

    steric_factor = 1
    alignment_factor = 1

    prec_settings = {  # settings for precomplex generator
        # Defaults:
        "num_spher": 2,
        "base_rot": 30,  #[next:60, rest:120]
        "planar_fix": True,
        "dist_weight": 1.0,
        "sadd_weight": 1.0,  # second add
        "max_min_weight": steric_factor * 1.0,
        "avg_dist_weight": steric_factor * 1.0,
        "add_dist_weight": alignment_factor * 1.0,
        "angle_diff_weight": alignment_factor * 1.0,
        "dummy_weight": alignment_factor * 6.0,
        "cylin_weight": steric_factor * 1.0,
        "numrot_weight": 1.0,  # not in use at the moment?
    }

    shuttle_mol = None
    shuttle=("None", 0, 0)
    solvent_pmg = None
    solvent_name = None


    with open(args.r, "r") as inp_data:
        if args.r[-4:] == "json":
            rule = json.load(inp_data)
        elif args.r[-3:] == "txt":
            rule = {"ADD": [], "BREAK":[]}
            lines = inp_data.readlines()
            for line in lines:
                if line.split()[0].lower() == "add":
                    if educt_names.index(line.split()[1]) == 1:
                        atom_id_1 = int(line.split()[2]) + natoms
                    else:
                        atom_id_1 = int(line.split()[2])
                    if educt_names.index(line.split()[3]) == 1:
                        atom_id_2 = int(line.split()[4]) + natoms
                    else:
                        atom_id_2 = int(line.split()[4])
                    rule["ADD"].append([educt_names.index(line.split()[1]), educt_names.index(line.split()[3]), atom_id_1, atom_id_2])
                elif line.split()[0].lower() == "break":
                    if educt_names.index(line.split()[1]) == 1:
                        atom_id_1 = int(line.split()[2]) + natoms
                        atom_id_2 = int(line.split()[3]) + natoms
                    else:
                        atom_id_1 = int(line.split()[2])
                        atom_id_2 = int(line.split()[3])
                    rule["BREAK"].append([educt_names.index(line.split()[1]), atom_id_1, atom_id_2])
                elif line.split()[0].lower() == "solvent":
                    solvent_name = line.split()[1]
                    shuttle = (solvent_name, line.split()[2], line.split()[3])
                else:
                    print(line.split()[0].lower())
                    sys.exit("ERROR in rule.txt file!")

    if solvent_name != None:
        smol = [{"name": solvent_name, "type": "xyz"}]
        shuttle_mol = getEducts(smol, verbose=True)
        solvent_pmg = Molecule.from_file(shuttle[0]+'.xyz')

    # Analyse reactants, define reaction vectors, etc.
    topoData = prep_topo(educt_names=educt_names, educts=educts, shuttle=shuttle, shuttle_mol=shuttle_mol)

    # Initialize preCompBuild class with the given educt structure, names, topology data and precomplex settings
    build_class = preCompBuild(educts=educts, educt_names=educt_names, structures=educts, topos=topoData,
                                prec_settings=prec_settings, solvent_pmg=solvent_pmg)

    # Generate Precomplexes
    precomplexes = generatePreComplex(rule=rule, build_class=build_class, educt_names=educt_names, educts=educts,
                                      shuttle=shuttle, shuttle_mol=shuttle_mol)

if __name__ == "__main__":
    main()
