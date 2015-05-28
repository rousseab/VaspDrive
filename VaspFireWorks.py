"""
 Implementation of FireWorks workflows for my own calculations.
"""

from fireworks.user_objects.firetasks.script_task import ScriptTask

from fireworks.core.firework import FireTaskBase, FWAction


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

        self.structure = d['structure']
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

