#================================================================================
# Suite of functions to analyse my VASP generated data
#================================================================================

import numpy as N 
import os
import sys
import json

import pymatgen

from pymatgen.serializers.json_coders import pmg_load
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.entries.compatibility import MaterialsProjectCompatibility

from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer


def compute_relative_energies(list_x_alkali,list_energies_per_unit):

    tol = 1e-8

    x_max = N.max(list_x_alkali)

    I0 = N.where( N.abs(list_x_alkali) < tol )
    E0 = N.min(list_energies_per_unit[I0])

    I1 = N.where( N.abs(list_x_alkali-x_max) < tol )[0]
    E1 = N.min(list_energies_per_unit[I1])

    list_relative_E = list_energies_per_unit-list_x/x_max*E1-(1.-list_x/x_max)*E0

    return list_relative_E 

def compute_voltage(list_x_alkali,list_energies_per_unit,E_alkali):
    """
    This function computes voltage plateaus for compound 
    """

    #  Build the lowest energy curve
    list_xV = []
    list_EV = []
    for x in N.sort(list(set(list_x_alkali))):
        I = N.where(list_x_alkali== x)[0]

        list_EV.append(N.min(list_energies_per_unit[I]))
        list_xV.append(x)

    list_EV = N.array(list_EV)
    list_xV = N.array(list_xV)

    list_V = E_alkali - (list_EV[1:]-list_EV[:-1])/(list_xV[1:]-list_xV[:-1])

    list_x_plateau = []
    list_V_plateau = []

    for i,V in enumerate(list_V):
        x0 = list_xV[i]
        x1 = list_xV[i+1]

        list_x_plateau.append(x0)
        list_x_plateau.append(x1)

        list_V_plateau.append(V)
        list_V_plateau.append(V)

    return N.array(list_x_plateau), N.array(list_V_plateau)


class AnalyseJsonData():
    """
    Class which will easily extract data from "run_data.json" files.
    The expected json files are in my own format, wrapping up most  of
    the relevant information in a VASP vasprun.xml and OUTCAR files. 
    """

    def __init__(self,list_json_data_filenames):

        self.list_json_data_filenames = list_json_data_filenames 

        self.parse_json_data()

        return

    def parse_json_data(self):
        self.list_data_dictionaries = []
        self.list_structures        = []


        for json_data_filename in self.list_json_data_filenames:
            with open(json_data_filename ,'r') as f:
                data_dictionary = json.load(f)
                self.list_data_dictionaries.append(data_dictionary) 

                structure_dict = data_dictionary['relaxation'][-1]['structure']
                structure = pymatgen.Structure.from_dict(structure_dict)

                self.list_structures.append(structure)

        return


    def extract_composition(self):

        list_compositions = []        
        for structure in self.list_structures:
            list_compositions.append(structure.composition)

        return list_compositions


    def extract_energies(self):

        list_energies = []
        for data_dictionary in self.list_data_dictionaries:
            energy = data_dictionary['relaxation'][-1]['electronic']['e_0_energy']
            list_energies.append(energy)

        return N.array(list_energies)



    def extract_magnetization(self,Element):
        MAG = []

        for data_dictionary, structure in zip(self.list_data_dictionaries,self.list_structures):
            list_d = data_dictionary['OUTCAR']['magnetization']                                                              
            list_mag = []

            for site, d in zip(structure,list_d):
                if site.specie == Element:
                    list_mag.append(d['tot'])

            MAG.append(list_mag)

        return MAG

class AnalyseMaterialsProjectJsonData():
    """
    Class which will wrap around boilerplate analysis of MaterialsProject-like
    json files, containing data extracted using borgs and queens.

    """

    def __init__(self):
        # some MP analysis power tools
        self.compat  = MaterialsProjectCompatibility()

        return

    def extract_alkali_energy(self, MP_alkali_json_data_filename):
        computed_entry = self._extract_MP_data(MP_alkali_json_data_filename)[0]
        processed_Alkali_entry = self.compat.process_entry(computed_entry)
        self.E_Alkali = processed_Alkali_entry.energy

        return

    def extract_phase_diagram_info(self,MP_phase_diagram_json_data_filename):

        computed_entries  = self._extract_MP_data(MP_phase_diagram_json_data_filename)
        processed_entries = self.compat.process_entries(computed_entries)

        pd = PhaseDiagram(processed_entries)
        self.phase_diagram_analyser = PDAnalyzer(pd)

        return

    def extract_processed_entries(self,MP_json_data_filename):

        computed_entries  = self._extract_MP_data(MP_json_data_filename)
        processed_entries = self.compat.process_entries(computed_entries)

        return processed_entries


    def extract_energies(self,MP_json_data_filename,alkali):

        processed_entries = self.extract_processed_entries(MP_json_data_filename)

        list_energy         = []
        list_alkali_content = []
        for entry in processed_entries:
            list_energy.append(entry.energy)
            list_alkali_content.append(entry.composition[alkali])

        list_energy         = N.array(list_energy)
        list_alkali_content = N.array(list_alkali_content )

        I = N.argsort(list_alkali_content )
        
        return list_alkali_content[I], list_energy[I]

    def _extract_MP_data(self,MP_data_filename):

        drone = VaspToComputedEntryDrone()
        queen = BorgQueen(drone, "dummy", 1)

        queen.load_data(MP_data_filename)
        computed_entries = queen.get_data()

        del drone
        del queen

        return computed_entries 
