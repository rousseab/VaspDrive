"""
 module to extract results from VASP computations
"""
import pymatgen
import os
import time
import shutil

from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen

from pymatgen.io.smartio import read_structure, write_structure
from pymatgen.io.vaspio import Outcar, Vasprun
from pymatgen.serializers.json_coders import pmg_dump


from pymatgen.io.vaspio_set import MPVaspInputSet
from pymatgen.io.vaspio.vasp_input import Incar, Potcar, Poscar
from pymatgen.entries.computed_entries import ComputedEntry

import ase.calculators.vasp as ase_vasp

Ha_to_eV = pymatgen.Ha_to_eV

def extract_DOS_data_with_ASE():
    """
    Routine tries to read DOSCAR file using ASE methods, then dumps 
    relevant results to json file.
    """

    data_dict = {}        
    try:
        ase_doscar = ase_vasp.VaspDos('DOSCAR')

        data_dict['energy']         = ase_doscar.energy
        data_dict['total_dos']      = ase_doscar.dos
        data_dict['integrated_dos'] = ase_doscar.integrated_dos


        list_projection = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
        list_spin       = ['up','down']

        good = True

        i_site = -1

        site_dos = []
        while good:
            i_site += 1

            site_dict = {}
            j_orbital = -1
            for p in list_projection:
                for s in list_spin:
                    j_orbital = +1
                    try:
                        dos_array = ase_doscar.site_dos(i_site,j_orbital)
                        site_dict[p+'-'+s] = dos_array 
                    except:
                        good = False

            if good:
                site_dos.append(site_dict)

        data_dict['site_dos']  = site_dos

        found_doscar = True
    except:
        print('DOSCAR file missing or not readable')
        found_doscar = False

    if found_doscar:
        pmg_dump(data_dict, 'ASE_DOS_data.json')

    return

def extract_json_data():
    """
    Routine tries to read VASP data into pymatgen objects, and
    then extracts only the relevant data. This is then written to json,
    allowing the voluminous OUTCAR and vasprun.xml files to be discarded.
    """

    try:
        o  = Outcar('OUTCAR')
        found_outcar = True
    except:
        print('OUTCAR file missing or not readable')
        found_outcar = False

    try:
        vr = Vasprun('vasprun.xml')
        found_vasprun = True
    except:
        print('vasprun.xml file missing or not readable')
        found_vasprun = False


    dictionary_data  = {}

    if found_outcar:
        dictionary_data['OUTCAR'] =  o.as_dict()

    if found_vasprun:
    
        try:
            # try to extract a Computed Entry object, using 
            # pymatgen technology
            drone = VaspToComputedEntryDrone()
            queen = BorgQueen(drone, './', 1)
            entry = queen.get_data()[0]

            dictionary_data['ComputedEntry'] =  entry.as_dict()
        except:
            print('ComputedEntry COULD NOT BE EXTRACTED BY PYMATGEN...')


        try:
            dictionary_data['DOS'] = vr.complete_dos.as_dict()
            pymatgen_dos_success = True
        except:
            print('DOS COULD NOT BE EXTRACTED BY PYMATGEN...')
            pymatgen_dos_success = False


        relaxation_data = []
        for step in vr.ionic_steps:

            data_dict = {}
            for key in ['forces','structure','stress']:
                if key in step:
                    data_dict[key] = step[key]

            data_dict['electronic'] = step['electronic_steps'][-1]

            relaxation_data.append(data_dict)

        dictionary_data['relaxation'] = relaxation_data 


    if found_outcar or found_vasprun:
        pmg_dump(dictionary_data, 'run_data.json')

    return    

def get_qstat_ids():
    os.system('qstat > qstat.log')
    with open('qstat.log','r') as f:
        lines = f.readlines()

    qstat_ids = map( lambda s: int(s.split()[0]), lines[2:])
    return qstat_ids 

class RepairJsonData(object):
    """ This simple class sets out to repair the run_data.json files
        which were generated without a ComputedEntry field, followed by 
        deletion of the vasprun.xml file.

        This class is unlikely to be very useful once I have repaired my 
        current pool of files, but let's keep it here for posterity.
    """

    def __init__(self,save_old_file=False):

        # read in what is there
        json_filename = 'run_data.json'

        try:
            with open(json_filename,'r') as f:
                data =json.load(f)
        except:
            print("CANNOT FIND run_data.json: CLASS MUST BE CALLED IN DIRECTORY WHERE FILE IS")
            return

        last_step = data['relaxation'][-1]

        structure = pymatgen.Structure.from_dict(last_step['structure'])
        composition = structure.composition
        energy = last_step['electronic']['e_0_energy']

        set = MPVaspInputSet()
        incar = set.get_incar(structure)
        poscar = set.get_poscar(structure)
        potcar = set.get_potcar(structure)

        parameters = self.parse_pseudo_inputs(poscar,incar,potcar)

        entry = pymatgen.entries.computed_entries.ComputedEntry(composition, energy, parameters=parameters )

        data['ComputedEntry'] = entry.to_dict()

        if save_old_file:
            shutil.move('run_data.json','run_data_OLD.json')

        pmg_dump(data, 'run_data.json')

        return


    def parse_pseudo_inputs(self,poscar,incar,potcar):
        """ 
        This code is copy-pasted from pymatgen's hive.py routine.
        It aims to generate run parameters from the known input, even
        if output was deleted. 
        """        
        param = {"hubbards":{}}
        if (poscar is not None and incar is not None and "LDAUU" in incar):
            param["hubbards"] = dict(zip(poscar.site_symbols, incar["LDAUU"]))

        param["is_hubbard"] = (
            incar.get("LDAU", False) and sum(param["hubbards"].values()) > 0
        ) if incar is not None else False
        param["run_type"] = None

        if incar is not None:
            param["run_type"] = "GGA+U" if param["is_hubbard"] else "GGA"

        potcar_symbols = []
        for p in potcar:
            potcar_symbols.append(p.header) 

        param["potcar_symbols"] = potcar_symbols

        return param
