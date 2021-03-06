"""
 Implementation of FireWorks workflows for my own calculations.
"""
import os
import shutil
import json
from copy import deepcopy
from fireworks.user_objects.firetasks.script_task import ScriptTask
from fireworks import Firework

from pymatgen import Structure

from fireworks.core.firework import FireTaskBase, FWAction
from VaspDrive.VaspSubmission import *
from VaspDrive.VaspExtract import *

from VaspDrive.VaspUStrategy import *

import re
import os
import sys

job_template = """#!/bin/bash
module purge
module add intel/2013u0
module add openmpi/1.6.3

VASP=/gpfsFS1/scratch/nobackup/projets/gc029/opt/vasp/4.6/vtstNC_XC_2013/vasp

source ~/.bash_profile
source /gpfsFS1/scratch/nobackup/projets/gc029/xw2738/BASH/create_machinefile.sh
create_machinefile
unset SGE_ROOT

time mpiexec -np $NBCORE -x LD_LIBRARY_PATH -machinefile machines $VASP &> output.txt
"""

BADER_job_template = """
module purge
module add intel/2013u0
module add openmpi/1.6.3

source ~/.bash_profile
source /gpfsFS1/scratch/nobackup/projets/gc029/xw2738/BASH/create_machinefile.sh
create_machinefile
unset SGE_ROOT

CHGSUM=/gpfsFS1/scratch/nobackup/projets/gc029/xw2738/TOOLS/PERL/chgsum.pl
BADER=/gpfsFS1/scratch/nobackup/projets/gc029/xw2738/depository/bader/bader

$CHGSUM AECCAR0 AECCAR2
$BADER CHGCAR -ref CHGCAR_sum
"""


class MyTestTask(FireTaskBase):
    """ Simple test task to play with Fireworks' various functionalities.
        Not directly related to Vasp per se.
    """

    _fw_name = "A Simple Test Task"

    def run_task(self, fw_spec):

        launch_dir = fw_spec['_launch_dir']
        counter    = fw_spec['counter']

        x = 0.05*counter
        with open('test_file_%i.txt'%counter,'w') as f:
            print >> f, 'x = %8.4f'%x

        if x < 0.5:
            new_counter = counter+1
            new_launch_dir = launch_dir.replace('relax_V%i'%counter,'relax_V%i'%new_counter)
            new_fw_spec = dict( _launch_dir = new_launch_dir,
                                counter = new_counter)

            new_fw = Firework(MyTestTask(), new_fw_spec )
            return FWAction(stored_data={'x': x}, additions=new_fw)

        else:            
            return FWAction()


    def find_iteration_number(self,launch_dir):
        
        last_string = launch_dir.split('/')[-1]

        number = int(re.search(r'\d+', last_string).group())

        return number

class MyVaspFireTask(FireTaskBase):
    """
    This task will write/hack the VASP inputs, so that they are ready to 
    roll.
    """

    _fw_name = 'MyVaspFireTask'

    def run_task(self, fw_spec):

        launch_dir = fw_spec['_launch_dir']

        self._load_params(fw_spec)

        if 'previous_launch_dir' in fw_spec:
            src_chg = fw_spec['previous_launch_dir']+'/CHGCAR'
            dst_chg = launch_dir+'/CHGCAR'

            src_wvf = fw_spec['previous_launch_dir']+'/WAVECAR'
            dst_wvf = launch_dir+'/WAVECAR'


        if self.job_type == 'coarse_relax' or self.job_type == 'relax' or \
                self.job_type == 'ground_state' or self.job_type == 'ramp':
            try:
                shutil.move(src_chg, dst_chg)
            except:
                print('CHGCAR could not be moved to working directory')
            try:
                shutil.move(src_wvf, dst_wvf)
            except:
                print('WAVECAR could not be moved to working directory')

            check = self.generate_VASP_inputs(self.structure, self.name, self.nproc, 
                                U_strategy_instance = self.U_strategy, 
                                supplementary_incar_dict = self.supplementary_incar_dict)

        elif self.job_type == 'DOS':
            try:
                shutil.copy(src_chg, dst_chg)            
            except:
                print('CHGCAR could not be COPIED to working directory')

            previous_vasp_dir = fw_spec['previous_launch_dir']
            self.get_MaterialsProject_DOS_VASP_inputs(self.structure, previous_vasp_dir, 
                                    self.nproc, kpoints_density=1000, 
                                    U_strategy_instance = self.U_strategy,
                                    supplementary_incar_dict = self.supplementary_incar_dict)

        # execute the BASH script. os.system should wait for the child process to end.
        os.system('bash job.sh')           

        # push an analysis job in the launchpad!
        new_fw = Firework(MyAnalysisFireTask(), fw_spec)

        # CASIR workaround
        self.change_gid()

        return FWAction(additions=new_fw)


    def change_gid(self):
        """ Casir is set up such that computation data should belong to group gc029, not to my user's group!
            This regularly crashes my jobs; here I will change ids every time an analysis takes place """

        uid_xw2738 = 1090
        gid_gc029 = 2005

        for dirpath, dirnames, filenames in os.walk('./'):
            for filename in filenames:
                path = os.path.join(dirpath,filename)
                os.chown(path,uid_xw2738,gid_gc029)

    def _load_params(self, d):
        """
        Extract parameters from dictionary, imposing default value for optional parameters
        if they are absent.
        """

        # It appears the structure gets serialized into a dictionary by the action of writing
        # the specs. Let's give it life again!
        self.structure = Structure.from_dict(d['structure'])

        self.name = d['name']
        self.job_type = d['job_type']
        self.version  = d['version']

        if 'nproc' in d:
            self.nproc = int(d['nproc'])
        else:
            self.nproc = 16

        if 'maximum_relaxations' in d:
            self.maximum_relaxations = int(d['maximum_relaxations'])
        else:
            self.maximum_relaxations = 9

        if 'distress_number_relaxations' in d:
            self.distress_number_relaxations = int(d['distress_number_relaxations'])
        else:
            self.distress_number_relaxations = 6


        if 'supplementary_incar_dict' in d:
            self.supplementary_incar_dict = d['supplementary_incar_dict']
        else:
            self.supplementary_incar_dict = None

        if 'strategy_type' in d:
            if  d['strategy_type'] == 'HexaCyanoFerrate':
                self.U_strategy = U_Strategy_HexaCyanoFerrate()

            elif  d['strategy_type'] == 'RAMP':
                if 'current_U_Fe' not in d:
                    print("==== ERROR: the variable 'current_U_Fe' must be defined for RAMPING!")
                    sys.exit()
                U_Fe = d['current_U_Fe']
                self.U_strategy = U_Strategy_RAMP(new_U_dict={'Fe':U_Fe})

            elif  d['strategy_type'] == 'MaterialsProject':
                self.U_strategy = U_Strategy_MaterialsProject(variable_magnetization_dict={'Fe':[5,4]})

            elif  d['strategy_type'] == 'MaterialsProject_V2':
                if 'variable_magnetization_dict' not in d:
                    print("==== ERROR: the variable 'variable_magnetization_dict' must be defined!")
                    sys.exit()
                variable_magnetization_dict = d['variable_magnetization_dict']
                self.U_strategy = U_Strategy_MaterialsProject_V2(variable_magnetization_dict)
            elif  d['strategy_type'] == 'U_Strategy_Yamada_Nitrogen':
                if 'variable_magnetization_dict' not in d:
                    print("==== ERROR: the variable 'variable_magnetization_dict' must be defined!")
                    sys.exit()
                variable_magnetization_dict = d['variable_magnetization_dict']
                self.U_strategy = U_Strategy_Yamada_Nitrogen(variable_magnetization_dict)
            else:
                print("UNKNOWN STRATEGY! FAIL HARD")
                sys.exit()
        else:
            self.U_strategy = None

    def generate_VASP_inputs(self, structure, name, nproc=16, U_strategy_instance = None, 
                                                            supplementary_incar_dict = None):
        """
        Inputs will be inspired by MaterialsProject, but this function is appropriate
        when we are seriously modifying the inputs such that they no longer conform to Materials Project.
        """

        # let's start with MP
        input_set = MPVaspInputSet()
        incar  = input_set.get_incar(structure)

        poscar_need_hack = False
        potcar_need_hack = False

        if U_strategy_instance != None:
            #  reading the structure here insures consistency, rather
            #  than having the strategy read the structure outside this driver.
            U_strategy_instance.read_structure(structure)

            # Generate all LDAU-related variables according to specified strategy.
            LDAU_dict, poscar_need_hack, potcar_need_hack = U_strategy_instance.get_LDAU()

            incar.update(LDAU_dict) 

        # set the number of parallel processors to sqrt(nproc), 
        # as recommended in manual.
        incar.update({'NPAR':int(np.sqrt(nproc))}) 

        if supplementary_incar_dict != None:
            incar.update(supplementary_incar_dict) 

        poscar = input_set.get_poscar(structure)
        kpoints = input_set.get_kpoints(structure)
        potcar = input_set.get_potcar(structure)

        incar.write_file('INCAR')
        poscar.write_file('POSCAR', vasp4_compatible = True)
        kpoints.write_file('KPOINTS')
        potcar.write_file('POTCAR')

        if poscar_need_hack:
            # do we need specialized hacking of the poscar because of the U strategy?       
            # For instance, if two Transition Metal sites with the same element 
            # must be treated as independent, Pymatgen is not capable of doing this
            # and we must HACK.
            new_poscar_lines = U_strategy_instance.get_new_poscar_lines()
            with open('POSCAR','w') as f:
                for line in new_poscar_lines:
                    print >>f, line.strip()

        if potcar_need_hack:
            # do we need specialized hacking of the potcar because of the U strategy?       
            new_potcar_symbols = U_strategy_instance.get_new_potcar_symbols(potcar)
            new_potcar = Potcar(new_potcar_symbols) 
            # overwrite the previous potcar
            new_potcar.write_file('POTCAR')

        with open('job.sh','w') as f:
            f.write(job_template)

        with open('clean.sh','w') as f:
            f.write(clean_template)

        return 0

    def get_MaterialsProject_DOS_VASP_inputs(self, structure, previous_vasp_dir, 
                                    nproc=16, kpoints_density=1000, U_strategy_instance = None, 
                                    supplementary_incar_dict=None):


        if supplementary_incar_dict != None:
            user_incar_dict = deepcopy(supplementary_incar_dict)
        else:
            user_incar_dict = {}


        poscar_need_hack = False
        potcar_need_hack = False
        if U_strategy_instance != None:
            #  reading the structure here insures consistency, rather
            #  than having the strategy read the structure outside this driver.
            U_strategy_instance.read_structure(structure)

            # Generate all LDAU-related variables according to specified strategy.
            LDAU_dict, poscar_need_hack, potcar_need_hack = U_strategy_instance.get_LDAU()

            user_incar_dict.update(LDAU_dict) 

        # set the number of parallel processors to sqrt(nproc), 
        # as recommended in manual.
        user_incar_dict.update({'NPAR':int(np.sqrt(nproc))}) 

        input_set = TetrahedronDosSet.from_previous_vasp_run(previous_vasp_dir,
                    kpoints_density=kpoints_density, user_incar_settings=user_incar_dict)

        if poscar_need_hack:
            # do we need specialized hacking of the poscar because of the U strategy?       
            # For instance, if two Transition Metal sites with the same element 
            # must be treated as independent, Pymatgen is not capable of doing this
            # and we must HACK.
            new_poscar_lines = U_strategy_instance.get_new_poscar_lines()
            with open('POSCAR','w') as f:
                for line in new_poscar_lines:
                    print >>f, line.strip()

        if potcar_need_hack:
            # do we need specialized hacking of the potcar because of the U strategy?       
            new_potcar_symbols = U_strategy_instance.get_new_potcar_symbols(potcar)
            new_potcar = Potcar(new_potcar_symbols) 
            # overwrite the previous potcar
            new_potcar.write_file('POTCAR')

        with open('job.sh','w') as f:
            f.write(job_template)

        with open('clean.sh','w') as f:
            f.write(clean_template)

        return 0

class MyAnalysisFireTask(FireTaskBase):
    """
    This task will read the VASP output, compute various metrics and decide what to do next.
    """

    _fw_name = 'MyAnalysisFireTask'

    def run_task(self, fw_spec):
        self.name = fw_spec['name']
        self.job_type = fw_spec['job_type']
        self.version  = fw_spec['version']

        # read the vasprun.xml and OUTCAR file,         
        # outputing the highlights to run_data.json
        extract_json_data()

        # get the structure that was written, and read it in
        structure, max_force = self.extract_run_data_info()

        # write a cif for posterity
        self.write_structure(structure)

        # Decide what to do next!
        firework_action = self.update_spec_and_launch_fireworks(structure, max_force, fw_spec)

        # change files' group id, to avoid CASIR idiosyncratic behavior.
        self.change_gid()

        return firework_action 

    def change_gid(self):
        """ Casir is set up such that computation data should belong to group gc029, not to my user's group!
            This regularly crashes my jobs; here I will change ids every time an analysis takes place """

        uid_xw2738 = 1090
        gid_gc029 = 2005

        for dirpath, dirnames, filenames in os.walk('./'):
            for filename in filenames:
                path = os.path.join(dirpath,filename)
                os.chown(path,uid_xw2738,gid_gc029)

            

    def define_job_dictionaries(self):
        """ Simple convenience to define job parameters for various situations """

        self.GS_dict  = dict(    EDIFF   =   1E-5,       # criterion to stop SCF loop, in eV
                                 ENCUT   =   520,        # HIGH encut
                                 PREC    =   'ACCURATE', # level of precision
                                 NSW     =     0,        # no ionic steps: fixed ions
                                 ICHARG  =     1,        # read in the CHGCAR file
                                 LORBIT  =   11,         # 11 prints out the DOS
                                 LCHARG  =   True,       # Write charge densities?
                                 LWAVE   =   False,      # write out the wavefunctions?
                                 NELM    =   100,        # maximum number of SCF cycles 
                                 ADDGRID =   True,       # fine FFT grid; not sure if this is needed for bader?
                                 LAECHG  =   True)       # Compute and write CORE electronic density, for BADER

        self.DOS_dict = dict(   EDIFF   =   1E-5,       # criterion to stop SCF loop, in eV
                                 ENCUT   =   520,       # HIGH encut
                                PREC    =   'ACCURATE', # level of precision
                                NSW     =   0,          # no ionic steps: fixed ions
                                LORBIT  =   11,         # 11 prints out the DOS
                                ADDGRID =  False,       # only necessary in ground state
                                LCHARG  =   False,      # Write charge densities?
                                LAECHG =    False,      # don't need all electron densities
                                LWAVE   =   False,      # write out the wavefunctions?
                                NELM    =   100,        # maximum number of SCF cycles 
                                ICHARG  =    11,        # non-scf density
                                IBRION  =    -1,        # don't move structure
                                ISMEAR  =    -5,        # tetrahedron integration
                                EMIN    =   -10,        # minimum energy for DOS
                                EMAX    =    10,        # maximum energy for DOS
                                NEDOS   =   2000)       # how many points for DOS calculation


        self.coarse_relax_dict = dict(    LCHARG  =  True,      # Write charge densities?
                                          PREC    = 'MEDIUM',   # level of precision
                                          EDIFF   =   1E-4,     # criterion to stop SCF loop, in eV
                                          EDIFFG  =  -5E-2,     # criterion to stop ionic relaxations. Negative means FORCES < |EDIFFG|
                                          ENCUT   =   260,      # LOW encut
                                          NELM    =    40,      # maximum number of SCF cycles 
                                          ICHARG  =     1,      # read in the CHGCAR file
                                          IBRION  =     2,      # use the CG algorithm
                                          ISIF    =     3,      # Relax cell shape 
                                          POTIM   =    0.5,    # controls step in relaxation algorithm
                                          NSW     =     30)     # max number of ionic steps: if it takes more, something is wrong.

        self.fine_relax_dict = dict(      LCHARG  =   True,       # Write charge densities?
                                          PREC    =   'ACCURATE', # level of precision
                                          EDIFF   =   1E-5,       # criterion to stop SCF loop, in eV
                                          EDIFFG  =  -5E-2,       # criterion to stop ionic relaxations. Negative means FORCES < |EDIFFG|
                                          ENCUT   =   520,        # HIGH encut
                                          NELM    =    40,        # maximum number of SCF cycles 
                                          ICHARG  =     1,        # read in the CHGCAR file
                                          IBRION  =     1,        # use the RMM-DIIS algorithm
                                          ISIF    =     3,        # Do relax cell shape 
                                          POTIM   =     0.5,      # controls step in relaxation algorithm
                                          NSW     =     30)       # max number of ionic steps: if it takes more, something is wrong.

        self.RAMP_dict = dict(   EDIFF   =   1E-5,       # criterion to stop SCF loop, in eV
                                 ENCUT   =   520,        # HIGH encut
                                 PREC    =   'ACCURATE', # level of precision
                                 NSW     =     0,        # no ionic steps: fixed ions
                                 ICHARG  =     1,        # read in the CHGCAR file
                                 LAECHG  =   False,      # DO NOT compute and write CORE electronic density
                                 LORBIT  =   11,         # 11 prints out the DOS
                                 LCHARG  =   True,       # Write charge densities?
                                 LWAVE   =   True,       # write out the wavefunctions?
                                 NELM    =   150,        # maximum number of SCF cycles 
                                 ISTART  =   1,          # read in wavefunctions, if they are present
                                 ALGO    =   'NORMAL')   # use block davidson, or else ZHEGV fails


    def update_spec_and_launch_fireworks(self,structure, max_force, fw_spec):

        fw_spec['structure'] = structure 
        fw_spec['previous_launch_dir'] = fw_spec['_launch_dir']


        if 'maximum_relaxations' in fw_spec:
            self.maximum_relaxations = int(fw_spec['maximum_relaxations'])
        else:
            self.maximum_relaxations = 9

        if 'distress_number_relaxations' in fw_spec:
            self.distress_number_relaxations = int(fw_spec['distress_number_relaxations'])
        else:
            self.distress_number_relaxations = 6


        formula  = structure.formula.replace(' ','')
        self.define_job_dictionaries()

        if self.job_type == 'relax': 
            if max_force < 0.05 or self.version >= self.maximum_relaxations:
                # relaxations are done, or it is hopless to reduce forces! Let's do a ground state!                    
                fw_spec['job_type'] = 'ground_state'
                fw_spec['name'] = formula+'_ground_state_V%i'%(self.version)  
                fw_spec['_launch_dir'] = fw_spec['top_dir']+'/ground_state_V%i/'%(self.version)  
                fw_spec['ground_state_dir'] = fw_spec['_launch_dir']  # for other post-processing to know where to get densities

                if 'supplementary_incar_dict' in fw_spec:
                    fw_spec['supplementary_incar_dict'].update(self.GS_dict) 
                else:
                    fw_spec['supplementary_incar_dict'] = self.GS_dict 

                new_fw = Firework(MyVaspFireTask(), fw_spec)

                return FWAction(additions=new_fw)

            else:
                #  we must relax some more
                fw_spec['job_type'] = 'relax'
                self.version += 1
                fw_spec['version'] = self.version
                fw_spec['name'] = formula+'_relax_V%i'%(self.version)  
                fw_spec['_launch_dir'] = fw_spec['top_dir']+'/relax_V%i/'%(self.version)  

                if self.version > self.distress_number_relaxations:
                    # relaxation is in distress!
                    potim = 0.2
                else:
                    # default value
                    potim = 0.5

                if 'supplementary_incar_dict' in fw_spec:
                    fw_spec['supplementary_incar_dict'].update(self.fine_relax_dict) 
                else:
                    fw_spec['supplementary_incar_dict'] = self.fine_relax_dict

                fw_spec['supplementary_incar_dict'].update({'potim':potim})

                new_fw = Firework(MyVaspFireTask(), fw_spec)
                return FWAction(additions=new_fw)

        if self.job_type == 'coarse_relax': 
            if max_force < 0.05 or self.version > 8: 
                # relaxations are done, or it is hopless to reduce forces! Let's go to fine_relax'
                fw_spec['job_type'] = 'relax'
                fw_spec['name'] = formula+'_relax_V%i'%(self.version)  
                fw_spec['_launch_dir'] = fw_spec['top_dir']+'/relax_V%i/'%(self.version)  
                fw_spec['ground_state_dir'] = fw_spec['_launch_dir']  # for other post-processing to know where to get densities

                if 'supplementary_incar_dict' in fw_spec:
                    fw_spec['supplementary_incar_dict'].update(self.fine_relax_dict)
                else:
                    fw_spec['supplementary_incar_dict'] = self.fine_relax_dict

                new_fw = Firework(MyVaspFireTask(), fw_spec)

                return FWAction(additions=new_fw)

            else:
                #  we must relax some more
                fw_spec['job_type'] = 'coarse_relax'
                self.version += 1
                fw_spec['version'] = self.version
                fw_spec['name'] = formula+'_coarse_relax_V%i'%(self.version)  
                fw_spec['_launch_dir'] = fw_spec['top_dir']+'/coarse_relax_V%i/'%(self.version)  

                if self.version > 6:
                    # relaxation is in distress!
                    potim = 0.2
                else:
                    # default value
                    potim = 0.5

                if 'supplementary_incar_dict' in fw_spec:
                    fw_spec['supplementary_incar_dict'].update(self.coarse_relax_dict) 
                else:
                    fw_spec['supplementary_incar_dict'] = self.coarse_relax_dict

                fw_spec['supplementary_incar_dict'].update({'potim':potim})

                new_fw = Firework(MyVaspFireTask(), fw_spec)
                return FWAction(additions=new_fw)

        if self.job_type == 'ramp': 
            # let's ramp up the U value!

            U_Fe = fw_spec['current_U_Fe']

            if U_Fe == 5.3:
                # we are done!
                return FWAction()

            self.version += 1
            fw_spec['version'] = self.version
            fw_spec['job_type'] = 'ramp'
            fw_spec['name'] = formula+'_ramp_V%i'%(self.version)  
            fw_spec['_launch_dir'] = fw_spec['top_dir']+'/ramp_V%i/'%(self.version)  

            if fw_spec['current_U_Fe'] < 0.1:
                #  go around the idiotic pymatgen behavior of reseting U = 5.3 if we set U = 0.
                fw_spec['current_U_Fe'] = 0.

            fw_spec['current_U_Fe'] += 0.1
            if 'supplementary_incar_dict' in fw_spec:
                fw_spec['supplementary_incar_dict'].update(self.RAMP_dict) 
            else:
                fw_spec['supplementary_incar_dict'] = self.RAMP_dict 

            new_fw = Firework(MyVaspFireTask(), fw_spec)
            return FWAction(additions=new_fw)

        if self.job_type == 'ground_state': 
            # let's do the DOS and BADER next!

            fw_spec_DOS = deepcopy(fw_spec)
            
            fw_spec_DOS['job_type'] = 'DOS'
            fw_spec_DOS['name'] = formula+'_DOS_V%i'%(self.version)  
            fw_spec_DOS['_launch_dir'] = fw_spec['top_dir']+'/DOS_V%i/'%(self.version)  


            if 'supplementary_incar_dict' in fw_spec_DOS:
                fw_spec_DOS['supplementary_incar_dict'].update(self.DOS_dict) 
            else:
                fw_spec_DOS['supplementary_incar_dict'] = self.DOS_dict

            new_DOS_fw = Firework(MyVaspFireTask(), fw_spec_DOS)

            fw_spec_BADER = deepcopy(fw_spec)
            fw_spec_BADER['job_type'] = 'Bader'
            fw_spec_BADER['name'] = formula+'_Bader_V%i'%(self.version)  
            fw_spec_BADER['_launch_dir'] = fw_spec['top_dir']+'/Bader_V%i/'%(self.version)  

            new_BADER_fw = Firework(MyBaderFireTask(), fw_spec_BADER)

            return FWAction(additions = [new_DOS_fw,new_BADER_fw] )

        if self.job_type == 'DOS' or self.job_type == 'Bader': 
            # we are done!
            return FWAction()

    def extract_run_data_info(self):
        """ read the json file and extract the structure """
        json_data_filename =  'run_data.json'
        file = open(json_data_filename ,'r')
        data_dictionary = json.load(file)
        file.close()

        last_step = data_dictionary['relaxation'][-1]
        structure_dict = last_step['structure']

        # we finally get the relaxed structure
        structure = Structure.from_dict(structure_dict)
        max_force = np.sqrt(np.sum(np.array(last_step['forces'])**2,axis=1)).max() 

        return structure, max_force 

    def write_structure(self,structure):

        formula  = structure.formula.replace(' ','')+'_%s_V%i'%(self.job_type,self.version)

        cif_filename = formula+'.cif'

        structure.to(fmt='cif',filename=cif_filename)

class MyBaderFireTask(FireTaskBase):
    """
    This task will move the AECAR files from the GS directory and compute the BADER charges.
    """

    _fw_name = 'MyBaderFireTask'

    def run_task(self, fw_spec):

        launch_dir = fw_spec['_launch_dir']
        gs_dir = fw_spec['ground_state_dir']

        # move the large density files to the workding directory
        for filename in ['AECCAR0','AECCAR2']:
            src = gs_dir+'/'+filename
            dst = launch_dir+'/'+filename
            shutil.move(src,dst)

        src = gs_dir+'/CHGCAR'
        dst = launch_dir+'/CHGCAR'
        shutil.copy(src,dst)

        with open('bader_job.sh','w') as f:
            f.write(BADER_job_template)

        os.system('bash bader_job.sh')           

        # CASIR workaround
        self.change_gid()

        # end of the line
        return FWAction()


    def change_gid(self):
        """ Casir is set up such that computation data should belong to group gc029, not to my user's group!
            This regularly crashes my jobs; here I will change ids every time an analysis takes place """

        uid_xw2738 = 1090
        gid_gc029 = 2005

        for dirpath, dirnames, filenames in os.walk('./'):
            for filename in filenames:
                path = os.path.join(dirpath,filename)
                os.chown(path,uid_xw2738,gid_gc029)
