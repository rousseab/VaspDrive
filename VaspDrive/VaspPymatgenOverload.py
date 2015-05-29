# coding: utf-8

from __future__ import division, unicode_literals

"""
This module defines the derived TetrahedronDosSet class, to facilitate the 
computation of DOS using the tetrahedron method. 
"""

__author__ = "Bruno Rousseau"
__copyright__ = "NONE"
__version__ = "NONE"
__maintainer__ = "Bruno Rousseau"
__email__ = "rousseau.bruno@gmail.com"
__date__ = "May 12, 2015"

from pymatgen.io.vaspio_set import MPNonSCFVaspInputSet, DictVaspInputSet, MODULE_DIR 

import os
import abc

import re
import traceback
import shutil
from functools import partial

import six
import numpy as np

from monty.serialization import loadfn

from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vaspio.vasp_output import Vasprun, Outcar
from pymatgen.serializers.json_coders import PMGSONable
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath


class TetrahedronDosSet(MPNonSCFVaspInputSet):
    """
    Implementation of VaspInputSet overriding MPNonSCFVaspInputSet,
    for non self-consistent field (NonSCF) calculation that follows
    a static run to calculate density of states (DOS) using the tetrahedron method.
    It is recommended to use the NonSCF from_previous_run method to construct
    the input set to inherit most of the functions.

    The prupose of this class is to overwrite the default contructor
    which prevents the use of the tetrahedron method in both the INCAR and the KPOINT
    files.
    Args:
        user_incar_settings (dict): A dict specify customized settings
            for INCAR. Must contain a NBANDS value, suggest to use
            1.2*(NBANDS from static run).
        mode: Line: Generate k-points along symmetry lines for
            bandstructure. Uniform: Generate uniform k-points
            grids for DOS.
        constrain_total_magmom (bool): Whether to constrain the total
            magmom (NUPDOWN in INCAR) to be the sum of the expected
            MAGMOM for all species. Defaults to False.
        kpoints_density (int): kpoints density for the reciprocal cell
            of structure. Might need to increase the default value when
            calculating metallic materials.
        kpoints_line_density (int): kpoints density to use in line-mode.
            Might need to increase the default value when calculating
            metallic materials.
        sort_structure (bool): Whether to sort structure. Defaults to
            False.
        sym_prec (float): Tolerance for symmetry finding
    """

    def __init__(self, user_incar_settings, 
                 constrain_total_magmom=False, sort_structure=False,
                 kpoints_density=1000, sym_prec=0.1):
        self.sym_prec = sym_prec
        
        DictVaspInputSet.__init__(
            self, "MaterialsProject Static",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")),
            constrain_total_magmom=constrain_total_magmom,
            sort_structure=sort_structure)

        self.user_incar_settings = user_incar_settings

        # impose tetrahedron method
        self.incar_settings.update(
            {"IBRION": -1, "ISMEAR": -5, "LCHARG": False,
             "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11})

        # this variable may have been used in the  ground state calculation,
        #  but is not relevant to the tetrahedron DOS calculation; its superfluous
        # presence in the INCAR leads to (human) confusion (VASP doesn't care).
        self.incar_settings.pop('SIGMA',None) 

        self.kpoints_settings.update({"kpoints_density": kpoints_density})

        # Set dense DOS output
        self.incar_settings.update({"NEDOS": 2001})

        if "NBANDS" not in user_incar_settings:
            raise KeyError("For NonSCF runs, NBANDS value from SC runs is "
                           "required!")
        else:
            self.incar_settings.update(user_incar_settings)

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               user_incar_settings=None,
                               copy_chgcar=True, make_dir_if_not_present=True,
                               kpoints_density=1000):
        """
        Generate a set of Vasp input files for NonSCF calculations from a
        directory of previous static Vasp run.

        Args:
            previous_vasp_dir (str): The directory contains the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            output_dir (str): The directory to write the VASP input files
                for the NonSCF calculations. Default to write in the current
                directory.
            user_incar_settings (dict): A dict specify customized settings
                for INCAR. 
            copy_chgcar (bool): Default to copy CHGCAR from SC run
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            kpoints_density (int): kpoints density for the reciprocal cell
                of structure. Might need to increase the default value when
                calculating metallic materials.
            kpoints_line_density (int): kpoints density to use in line-mode.
                Might need to increase the default value when calculating
                metallic materials.
        """
        user_incar_settings = user_incar_settings or {}

        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=None)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
            previous_incar = vasp_run.incar
        except:
            traceback.print_exc()
            raise RuntimeError("Can't get valid results from previous run. prev dir: {}".format(previous_vasp_dir))

        #Get a Magmom-decorated structure
        structure = TetrahedronDosSet.get_structure(vasp_run, outcar,
                                                       initial_structure=True)

        nscf_incar_settings = TetrahedronDosSet.get_incar_settings(vasp_run,
                                                                      outcar)

        mpnscfvip = TetrahedronDosSet(nscf_incar_settings, kpoints_density=kpoints_density)

        mpnscfvip.write_input(structure, output_dir, make_dir_if_not_present)

        if copy_chgcar:
            try:
                shutil.copyfile(os.path.join(previous_vasp_dir, "CHGCAR"),
                                os.path.join(output_dir, "CHGCAR"))
            except Exception as e:
                traceback.print_exc()
                raise RuntimeError("Can't copy CHGCAR from SC run" + '\n'
                                   + str(e))

        #Overwrite necessary INCAR parameters from previous runs
        #  this is already done in hte __init__; is it necessary here?
        previous_incar.update({"IBRION": -1, "ISMEAR": -5, 
                               "LCHARG": False, "LORBIT": 11, "LWAVE": False,
                               "NSW": 0, "ISYM": 0, "ICHARG": 11})

        previous_incar.update(nscf_incar_settings)
        previous_incar.update(user_incar_settings)

        previous_incar.pop("MAGMOM", None)
        previous_incar.pop('SIGMA',None) 

        previous_incar.write_file(os.path.join(output_dir, "INCAR"))

        # Perform checking on INCAR parameters
        if any([previous_incar.get("NSW", 0) != 0,
                previous_incar["IBRION"] != -1,
                previous_incar["ICHARG"] != 11,
               any([sum(previous_incar["LDAUU"]) <= 0,
                    previous_incar["LMAXMIX"] < 4])
               if previous_incar.get("LDAU") else False]):
            raise ValueError("Incompatible INCAR parameters!")

    def get_kpoints(self, structure):
        """
        Get a KPOINTS file for NonSCF calculation.  kpoints are
        Gamma-centered mesh grid. 

        Args:
            structure (Structure/IStructure): structure to get Kpoints
        """
        kppa = self.kpoints_settings["kpoints_density"]
        kpoints = Kpoints.automatic_gamma_density(structure, kppa)

        return kpoints

    def write_input(self, structure, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        """
        Writes a set of VASP input to a directory.

        Must be overloaded to be compatible with VASP 4.6.

        Args:
            structure (Structure/IStructure): Structure to write VASP input
                files for.
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            include_cif (bool): Whether to write a CIF file in the output
                directory for easier opening by VESTA.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.get_all_vasp_input(structure).items():

            if k == 'POSCAR':
                v.write_file(os.path.join(output_dir, k), vasp4_compatible = True)

            else:
                v.write_file(os.path.join(output_dir, k))

            if k == "POSCAR" and include_cif:
                v.structure.to(
                    filename=os.path.join(output_dir,
                                          "%s.cif" % v.structure.formula))


