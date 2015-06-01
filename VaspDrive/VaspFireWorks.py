"""
 Implementation of FireWorks workflows for my own calculations.
"""

from fireworks.user_objects.firetasks.script_task import ScriptTask
from fireworks import Firework

from pymatgen import Structure

from fireworks.core.firework import FireTaskBase, FWAction
from VaspSubmission import *
import re

class BuildVaspInputTask(FireTaskBase):
    """ FireWork task which wraps around the generation of VASP inputs """

    required_params = ["structure", "workdir", "job_name"]
    optional_params = ["nproc", "supplementary_incar_dict"]

    _fw_name = 'BuildVaspInputTask'

    def run_task(self, fw_spec):
        """ overload the run_task method to write the Vasp input """
        self._load_params(self)

        check = get_MaterialsProject_VASP_inputs(self.structure, self.workdir, self.job_name, 
                                    nproc=self.nproc, supplementary_incar_dict = self.supplementary_incar_dict)

    def _load_params(self, d):
        """
        Extract parameters from dictionary, imposing default value for optional parameters
        if they are absent.
        """

        # It appears the structure gets serialized into a dictionary by the action of writing
        # the specs. Let's give it life again!
        self.structure = Structure.from_dict(d['structure'])

        self.workdir = d['workdir']
        self.job_name = d['job_name']

        if 'nproc' in d:
            self.nproc = int(d['nproc'])
        else:
            self.nproc = 16

        if 'supplementary_incar_dict' in d:
            self.supplementary_incar_dict = d['supplementary_incar_dict']
        else:
            self.supplementary_incar_dict = None



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
