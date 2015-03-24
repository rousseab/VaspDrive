#================================================================================
# Suite of functions to analyse my VASP generated data
#================================================================================

import numpy as N 
import os
import sys

import pymatgen
from pymatgen.serializers.json_coders import pmg_load

import json



class AnalyseJsonData():
    """
    Class which will easily extract data from "run_data.json" files.
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

