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
    
    def __init__(self):
        """ 
        Don't do anything; this object is meant to be passed to a driver, which 
        will decide what to do.
        """

        self._LDAU_KEYS = ['LDAUTYPE', 'LDAUPRINT', 'MAGMOM', 'LDAUL', 'LDAUJ', 'LDAUU', 'LDAU'] 


    def read_structure(self,structure):
        """ 
        Get a copy of the structure object.         
        """

        # let's make sure we don't modifiy the input structure.
        # Side effects are bad practice!
        self.structure = structure.copy()

    def get_LDAU(self):
        """ Produce LDAU related variables, to be passed to VASP as strings """

        # let's simply use the Materials Project results as default.

        dummy_input_set = MPVaspInputSet()
        dummy_incar  = dummy_input_set.get_incar(self.structure)

        LDAU_dict = {}
        
        for key in self._LDAU_KEYS: 
            if key in dummy_incar:
                LDAU_dict[key] = dummy_incar[key]

        # no need to hack  the poscar or potcar
        poscar_need_hack = False
        potcar_need_hack = False

        return LDAU_dict, poscar_need_hack, potcar_need_hack  

    def hack_poscar(self):
        """ nothing to do!"""
        return

    def hack_potcar(self):
        """ nothing to do!"""
        return


class U_Strategy_HexaCyanoFerrate(U_Strategy):
    """
    Derived Class to treat specifically the case of hexacyanoferrate, where
    we want to impose two different values of U on Fe, depending on its neighbors.
    """

    def get_LDAU(self, U_Fe_N = 7., U_Fe_C = 3.):
        """ Overload this method to deal with the specific HexaCyanoFerrate case. 
            Different U values will be used for different Fe environments.
        """

        Fe = pymatgen.Element('Fe')

        Fe_hi = pymatgen.Specie('Fe',oxidation_state=4)
        Fe_lo = pymatgen.Specie('Fe',oxidation_state=1)

        N = pymatgen.Element('N')
        C = pymatgen.Element('C')

        neibhor_distance = 2.6 # angstrom
        # Spoof pymatgen by substituting high spin iron and low spin iron.
        for i,site in enumerate(self.structure.sites):
            if site.specie == Fe:
                neighbor = self.structure.get_neighbors(site,neibhor_distance)[0][0]

                if neighbor.specie == N:
                    self.structure.replace(i,Fe_hi)
                elif neighbor.specie == C:
                    self.structure.replace(i,Fe_lo)

        # Sort structure, so that decorated sites are 
        # next to each other
        self.structure.sort()
        
        # Force GGA+U on Fe site
        LDAUJ = ''
        LDAUL = ''
        LDAUU = ''
        MAGMOM = ''


        species_dict = OrderedDict()
        for s in self.structure.types_of_specie:
            species_dict[s] = 0.

            LDAUJ += ' 0'
            if s == Fe_lo:
                LDAUL += ' 2'
                LDAUU += ' %2.1f'%U_Fe_C
            elif s == Fe_hi:
                LDAUL += ' 2'
                LDAUU += ' %2.1f'%U_Fe_N
            else:
                LDAUL += ' 0'
                LDAUU += ' 0'

        for s in self.structure.sites:
            species_dict[s.specie] += 1.

        for s in self.structure.types_of_specie:

            if s == Fe_lo:
                MAGMOM += ' %i*1'%species_dict[s]  # low spin
            elif s == Fe_hi:
                MAGMOM += ' %i*5'%species_dict[s]  # high spin

            else:
                MAGMOM += ' %i*0.6'%species_dict[s] 

        LDAU_dict = { 'LDAU':True       # use LDA+U (GGA+U in fact)
                      'LDAUTYPE':2,     # simplified Dudarev Formalism
                      'LDAUPRINT':1,    # talk to me
                      'MAGMOM':MAGMOM,  # magnetic moments
                      'LDAUL':LDAUL, 
                      'LDAUJ':LDAUJ, 
                      'LDAUU':LDAUU}

        poscar_need_hack = True
        potcar_need_hack = True

        return  LDAU_dict, poscar_need_hack, potcar_need_hack 
