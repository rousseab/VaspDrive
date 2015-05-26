"""
   Module defining the U_Strategy class, a flexible object to
   generate U related information in a way which does not rely
   on the ordering in the CIF file.
"""

import pymatgen
from pymatgen.io.vaspio_set import MPVaspInputSet

import numpy as np


class U_Strategy(object):
    """
    Class to create LDAU strings for VASP input, beyond what the MaterialsProject already
    provide.
    """
    
    def __init__(self,structure):

        # let's make sure we don't modifiy the input structure.
        # Side effects are bad practice!
        self.structure = structure.copy()

        self._LDAU_KEYS = ['LDAUTYPE', 'LDAUPRINT', 'MAGMOM', 'LDAUL', 'LDAUJ', 'LDAUU', 'LDAU'] 


    def get_LDAU_strings(self):
        """ Produce LDAU related variables, to be passed to VASP as strings """

        # let's simply use the Materials Project results as default.

        dummy_input_set = MPVaspInputSet()
        dummy_incar  = dummy_input_set.get_incar(self.structure)

        LDAU_dict = {}
        
        for key in self._LDAU_KEYS: 
            if key in dummy_incar:
                LDAU_dict[key] = dummy_incar[key]

        return LDAU_dict

