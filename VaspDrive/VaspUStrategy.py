"""
   Module defining the U_Strategy class, a flexible object to
   generate U related information in a way which does not rely
   on the ordering in the CIF file.
"""

from copy import deepcopy
import pymatgen
from collections import OrderedDict, namedtuple
from pymatgen.io.vaspio_set import MPVaspInputSet
from pymatgen.io.vaspio.vasp_input import Poscar, Potcar
import sys


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

        self.structure_has_been_read = False


    def read_structure(self,structure):
        """ 
        Get a copy of the structure object.         
        """
        # let's make sure we don't modifiy the input structure.
        # Side effects are bad practice!
        self.structure = structure.copy()
        self.structure_has_been_read = True


    def check_structure_is_read(self):
        if not self.structure_has_been_read: 
            print('NEED TO READ STRUCTURE BEFORE PROCEEDING FURTHER!')
            sys.exit()


    def get_LDAU(self):
        """ Produce LDAU related variables, to be passed to VASP as strings """

        # let's simply use the Materials Project results as default.

        self.check_structure_is_read()

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

    def get_new_poscar_lines(self):
        """ nothing to do!"""
        return

    def get_new_potcar_symbols(self,old_potcar):
        """ nothing to do!"""
        return


class U_Strategy_MaterialsProject(U_Strategy):
    """
    Class to create LDAU strings for VASP input, starting from what the MaterialsProject already
    provides, but with variable starting magnetization on the transition metal sites.

    Who to reduce first? This is a tricky problem, which would be difficult to solve robustly and generally.
    The following is good enough "right now", but may warrant improvement in the future.
    """
    
    def __init__(self,variable_magnetization_dict={'Fe':[5,4]}):
        """ 

        input:        
            variable_magnetization_dict: 
                                    - key: elements which can have different oxydation states
                                    - value: magnetization to be put in MAGMOM, most oxidized case first.

        Note: 
            In a perfect world, we would import a strong AI (or just a better code) which would automatically
            determine the best starting MAGMOM for a given partially sodiated structure. In this world, however,            
            I don't have time to design/build/test the perfect solution and some a priori info will have to be passed to this
            object.

        """

        self.structure_has_been_read = False

        self._LDAU_KEYS = ['LDAUTYPE', 'LDAUPRINT', 'MAGMOM', 'LDAUL', 'LDAUJ', 'LDAUU', 'LDAU'] 

        self.variable_magnetization_dict = variable_magnetization_dict 


    def sort_TM_sites_by_Na_distance(self,Na_indices):
        """
        Identify the next item to be reduced, with the rules:
            - reduce Fe-C before Fe-N
            - reduce Fe nearest Na first
        """

        list_oxidizable_site_indices = []        

        distance_table = self.structure.distance_matrix[:,Na_indices]


        # find all elements in the structure which can be oxidized,
        # and the minimum distance to a Na atom
        list_d = []
        for i_s, site in enumerate(self.structure.sites):
            if site.specie.symbol in self.variable_magnetization_dict: 
                dist = np.min(distance_table[i_s])

                list_oxidizable_site_indices.append(i_s)
                list_d.append(dist) 

        # Sort this list_according to distance
        I = np.argsort(list_d)
        list_oxidizable_site_indices = np.array(list_oxidizable_site_indices)[I]

        return list_oxidizable_site_indices 
  
    def build_magmom(self,list_oxidizable_site_indices,number_of_electrons):
        """
        Build MAGMOM, given that some sites must be reduced
        """

        MAGMOM = []
        ne = number_of_electrons

        dict_oxidizable = {}

        for i_s in list_oxidizable_site_indices:
            symbol = self.structure.sites[i_s].specie.symbol
            if ne > 0:
                dict_oxidizable[i_s] = self.variable_magnetization_dict[symbol][1]
                ne -= 1
            elif ne == 0:
                dict_oxidizable[i_s] = self.variable_magnetization_dict[symbol][0]

            else:
                print("SOMETHING IS WRONG. REVIEW CODE!")
                sys.exit()

        for i_s, site in enumerate(self.structure):
            if i_s in dict_oxidizable:
                MAGMOM.append(dict_oxidizable[i_s])
            else:
                MAGMOM.append(0.6)

        return MAGMOM 



    def get_LDAU(self):
        """ Produce LDAU related variables, to be passed to VASP as strings """

        # let's simply use the default as a first step

        LDAU_dict, poscar_need_hack, potcar_need_hack = super(U_Strategy_MaterialsProject, self).get_LDAU()

        Na_indices = self.structure.indices_from_symbol('Na')
        number_of_electrons = len(Na_indices) 

        if number_of_electrons == 0:
            # structure is fully oxidized, no need to do anything
            # no need to hack  the poscar or potcar
            poscar_need_hack = False

        else:
            # structure must be partially reduced; hack MAGMOM
            list_oxidizable_site_indices = self.sort_TM_sites_by_Na_distance(Na_indices)

            MAGMOM = self.build_magmom(list_oxidizable_site_indices,number_of_electrons)
            LDAU_dict['MAGMOM'] = MAGMOM 

            poscar_need_hack = True

        return LDAU_dict, poscar_need_hack, potcar_need_hack  


class U_Strategy_HexaCyanoFerrate(U_Strategy):
    """
    Derived Class to treat specifically the case of hexacyanoferrate, where
    we want to impose two different values of U on Fe, depending on its neighbors.
    """

    def __init__(self):

        # create a boolean "state" variable which
        # will be used to make sure the structure has been modified
        # to account for different Fe sites having different U.
        # This will serve to impose routines being called in the right order.
        self.structure_has_been_read = False
        self.structure_has_been_modified = False


    def check_structure_is_modified(self):
        """ I'm sure this can be done better with @property..."""
        if not self.structure_has_been_modified: 
            print('NEED TO MODIFY STRUCTURE BEFORE PROCEEDING FURTHER!')
            sys.exit()

    def _modify_structure(self):

        # We'll use "oxidation_state" as a stand-in for "magnetisation"
        # surrounded by nitrogen
        self.Fe_HS_2plus = pymatgen.Specie('Fe',oxidation_state=4)
        self.Fe_HS_3plus = pymatgen.Specie('Fe',oxidation_state=5)

        # surrounded by carbon
        self.Fe_LS_2plus = pymatgen.Specie('Fe',oxidation_state=0)
        self.Fe_LS_3plus = pymatgen.Specie('Fe',oxidation_state=1)

        Fe = pymatgen.Element('Fe')
        Na = pymatgen.Element('Na')
        N = pymatgen.Element('N')
        C = pymatgen.Element('C')


        Na_indices = self.structure.indices_from_symbol('Na')
        
        neibhor_distance = 2.6 # angstrom

        # Interrogate the structure: for each Fe, what is the site index, the nearest neighbor
        # and the distance to Na?

        # create a nice verbose object for this job
        Iron_Ion = namedtuple('Iron_Ion',['structure_index', 'neighbor', 'distance_to_Na'])

        # Count the number of reducing electrons which must be introduced
        number_of_reducing_electrons = self.structure.composition['Na']


        list_Fe_ions = []

        # Identify all the Fe ions in the structure 
        for i,site in enumerate(self.structure.sites):
            if site.specie == Fe:
                neighbor = self.structure.get_neighbors(site,neibhor_distance)[0][0].specie

                if number_of_reducing_electrons > 0:
                    distance = np.min(self.structure.distance_matrix[i,Na_indices])
                else:
                    distance = np.infty

                list_Fe_ions.append( Iron_Ion(i,neighbor,distance) )

        # Spoof pymatgen by substituting reduced iron for the plain vanilla Fe.
        while number_of_reducing_electrons > 0:
            next_ion = self.find_next_site_to_reduce_pop(list_Fe_ions)

            if next_ion.neighbor == N:
                self.structure.replace(next_ion.structure_index,self.Fe_HS_2plus )
            elif next_ion.neighbor == C:
                self.structure.replace(next_ion.structure_index,self.Fe_LS_2plus )

            number_of_reducing_electrons -= 1

        # Spoof pymatgen by substituting oxidized iron for the plain vanilla Fe for the remaining sites
        for next_ion in list_Fe_ions:
            if next_ion.neighbor == N:
                self.structure.replace(next_ion.structure_index,self.Fe_HS_3plus )
            elif next_ion.neighbor == C:
                self.structure.replace(next_ion.structure_index,self.Fe_LS_3plus )

        # Sort structure, so that decorated sites are 
        # next to each other
        self.structure.sort()
        self.structure_has_been_modified = True
        return


    def find_next_site_to_reduce_pop(self,list_Fe_ions):
        """
        Identify the next item to be reduced, with the rules:
            - reduce Fe-C before Fe-N
            - reduce Fe nearest Na first
        """
        N = pymatgen.Element('N')
        C = pymatgen.Element('C')

        # Are there Fe-C ions?
        sub_list_Fe_ions = []
        for ion in list_Fe_ions:
            if ion.neighbor == C:
                sub_list_Fe_ions.append(ion)

        if len(sub_list_Fe_ions) == 0:
            # there are no Fe-C ion. They must all be Fe-N
            sub_list_Fe_ions = deepcopy(list_Fe_ions)

        # Find the site with the smallest Na distance
        sub_list_Na_d = []
        for ion in sub_list_Fe_ions:
            sub_list_Na_d.append(ion.distance_to_Na)

        next_ion = sub_list_Fe_ions[ np.argmin(sub_list_Na_d) ]


        # identify index of the identified ion, and 
        # pop it out of the list
        for i, ion in enumerate(list_Fe_ions):
            if ion == next_ion:
                list_Fe_ions.pop(i)
                break

        return next_ion


    def get_LDAU(self, U_Fe_N = 7., U_Fe_C = 5.):
        """ Overload this method to deal with the specific HexaCyanoFerrate case. 
            Different U values will be used for different Fe environments.
        """

        self.check_structure_is_read()

        self._modify_structure()


        # Initialize various strings
        LDAUJ = ''
        LDAUL = ''
        LDAUU = ''
        MAGMOM = ''

        # count the number of each species. 
        #   NOTE: pymatgen.composition does not understand decorations;
        #         we have to hack this way.

        self.species_dict = OrderedDict()

        for s in self.structure.types_of_specie:
            self.species_dict[s] = 0.

        for s in self.structure.sites:
            self.species_dict[s.specie] += 1.

        # Generate MAGMOM, which must distinguish every magnetization state
        for s in self.structure.types_of_specie:

            if s == self.Fe_HS_2plus:
                MAGMOM += ' %i*4'%self.species_dict[s]  # low spin
            elif s == self.Fe_HS_3plus:
                MAGMOM += ' %i*5'%self.species_dict[s]  # high spin
            elif s == self.Fe_LS_2plus:
                MAGMOM += ' %i*0'%self.species_dict[s]  # low spin
            elif s == self.Fe_LS_3plus:
                MAGMOM += ' %i*1'%self.species_dict[s]  # high spin
            else:
                MAGMOM += ' %i*0.6'%self.species_dict[s] 


        for s in self.structure.types_of_specie:
            LDAUJ += ' 0'

            if s == self.Fe_LS_2plus or s == self.Fe_LS_3plus: 
                LDAUL += ' 2'
                LDAUU += ' %2.1f'%U_Fe_C
            elif s == self.Fe_HS_2plus or s == self.Fe_HS_3plus: 
                LDAUL += ' 2'
                LDAUU += ' %2.1f'%U_Fe_N
            else:
                LDAUL += ' 0'
                LDAUU += ' 0'


        LDAU_dict = { 'LDAU':True,      # use LDA+U (GGA+U in fact)
                      'LDAUTYPE':2,     # simplified Dudarev Formalism
                      'LDAUPRINT':1,    # talk to me
                      'MAGMOM':MAGMOM,  # magnetic moments
                      'LDAUL':LDAUL, 
                      'LDAUJ':LDAUJ, 
                      'LDAUU':LDAUU,
                      'LMAXMIX':4} # this is essential for d-element GGA+U, but not the VASP default

        poscar_need_hack = True
        potcar_need_hack = True

        return  LDAU_dict, poscar_need_hack, potcar_need_hack 

    def get_new_poscar_lines(self):
        """
        Pymatgen does not realize that decorated elements (species) should be
        treated as different sites in VASP.

        This routine hacks the poscar file (which has already been hacked to go from 5.x to 4.6)
        to include the correct number of distinct sites.
        """

        self.check_structure_is_modified()

        # get a poscar consistent with the internally modified structure
        poscar = Poscar(self.structure)

        lines = poscar.get_string(vasp4_compatible=True).split('\n')
        
        hack_line = ''
        for element, number in self.species_dict.items():
            hack_line += ' %i'%number

        lines[5] = hack_line

        return lines


    def get_new_potcar_symbols(self, old_potcar):
        """ Get potcar symbols, accounting for repetition of Fe in the structure"""
    
        self.check_structure_is_modified()

        new_potcar_symbols = []

        for element, number in self.species_dict.items():

            symbol = element.symbol
            # the element symbol will certainly be in the potcar symbol
            # find the closest condender and add to the new list;
            # this should create duplicates.

            # Careful! This is a bit tricky 'N' is in 'Na', so a better match
            # algorithm is necessary
            Found = False
            for psymb in old_potcar.symbols:
                if symbol == psymb.split('_')[0]:
                    new_potcar_symbols.append(psymb) 
                    Found = True
                    break
            if not Found:        
                # Throw a fit
                print('PSEUDOPOTENTIAL FOR %s NOT FOUND'%symbol )
                print('this most certainly indicates a weakness in the algorithm used;')
                print('review code; for now, FAIL HARD')
                sys.exit()

        return new_potcar_symbols


class U_Strategy_HexaCyanoFerrate_U_is_5_7(U_Strategy_HexaCyanoFerrate):
    """
    Derived Class to treat specifically the case of hexacyanoferrate, where
    we want to impose two different values of U on Fe, depending on its neighbors.

    This derived class simply changes the default values of the U parameters.
    """

    def get_LDAU(self):
        super(U_Strategy_HexaCyanoFerrate_U_is_5_7, self).get_LDAU(U_Fe_N = 7., U_Fe_C = 5.)
