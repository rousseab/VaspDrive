"""
 Implementation of FireWorks workflows for my own calculations.
"""

from fireworks.user_objects.firetasks.script_task import ScriptTask
from fireworks import Firework

from pymatgen import Structure

from fireworks.core.firework import FireTaskBase, FWAction
from VaspDrive.VaspSubmission import *

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

class MyTestTask(FireTaskBase):
    """ Simple test task to play with Fireworks' various functionalities.
        Not directly related to Vasp per se.
    """

    _fw_name = "A Simple Test Task"

    def run_task(self, fw_spec):

        launch_dir = fw_spec['_launch_dir']
        iteration_number = self.find_iteration_number(launch_dir.strip('/'))

        x = np.random.random(1)[0]
        with open('test_file_%i.txt'%iteration_number,'w') as f:
            print >> f, 'x = %8.4f'%x

        if x < 0.8:
            new_launch_dir = launch_dir.replace('%i'%iteration_number,'%i'%(iteration_number+1))
            new_fw_spec = dict( _launch_dir = new_launch_dir) 

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

        self.generate_VASP_inputs(self.structure, self.job_name, self.nproc, 
                U_strategy_instance = self.U_strategy, supplementary_incar_dict = self.supplementary_incar_dict)

        # execute the BASH script
        os.system('bash job.sh')           



    def _load_params(self, d):
        """
        Extract parameters from dictionary, imposing default value for optional parameters
        if they are absent.
        """

        # It appears the structure gets serialized into a dictionary by the action of writing
        # the specs. Let's give it life again!
        self.structure = Structure.from_dict(d['structure'])
        self.job_name = d['job_name']

        if 'nproc' in d:
            self.nproc = int(d['nproc'])
        else:
            self.nproc = 16

        if 'supplementary_incar_dict' in d:
            self.supplementary_incar_dict = d['supplementary_incar_dict']
        else:
            self.supplementary_incar_dict = None

        if 'strategy_type' in d:
            if  d['strategy_type'] == 'HexaCyanoFerrate':
                self.U_strategy = U_Strategy_HexaCyanoFerrate()
            else:
                print("UNKNOWN STRATEGY! FAIL HARD")
                sys.exit()
        else:
            self.U_strategy = None
        

    def generate_VASP_inputs(self, structure, job_name, nproc=16, U_strategy_instance = None, supplementary_incar_dict = None):
        """
        Inputs will be inspired by MaterialsProject, but this function is appropriate
        when we are seriously modifying the inputs such that they no longer conform to Materials Project.
        """

        # let's start with MP
        input_set = MPVaspInputSet()
       
        incar  = input_set.get_incar(structure)

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
