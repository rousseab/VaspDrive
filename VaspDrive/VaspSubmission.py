"""
   Module wrapping useful pymatgen stuff to run VASP computations.
"""

import pymatgen
import numpy as np
import os

from pymatgen import read_structure, write_structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Kpoints, Potcar,  PotcarSingle
from pymatgen.io.vaspio.vasp_output import Vasprun

from pymatgen.io.vaspio_set import MPVaspInputSet, MPNonSCFVaspInputSet

from pymatgen import Structure

from VaspPymatgenOverload import TetrahedronDosSet

Ha_to_eV = pymatgen.Ha_to_eV


PSEUDOPOTENTIALS_DIRECTORY = '/gpfsFS1/scratch/nobackup/projets/gc029/xw2738/VASP/PSEUDOPOTENTIALS/'

#===========================================
# Template for submission files
#===========================================
submit_template = """#!/bin/bash
#$ -N {0}
#$ -pe dist16ppn {1}
#$ -cwd
#$ -j y
#$ -m n
#$ -p -1023
#$ -S /bin/bash
# Commands before execution
module purge
module add intel/2013u0
module add openmpi/1.6.3

VASP=/gpfsFS1/scratch/nobackup/projets/gc029/opt/vasp/4.6/vtstNC_XC_2013/vasp

source ~/.bash_profile
source /gpfsFS1/scratch/nobackup/projets/gc029/xw2738/BASH/create_machinefile.sh
create_machinefile
unset SGE_ROOT

time mpiexec -np $NBCORE -x LD_LIBRARY_PATH -machinefile machines $VASP &> output.txt

# remove large files 
#rm -rf  OUTCAR output.txt
"""

clean_template = """
rm -rf CHG DOSCAR IBZKPT OSZICAR output.txt res_out.txt vasprun.xml XDATCAR CHGCAR CONTCAR EIGENVAL machines OUTCAR PCDAT WAVECAR
"""


#===========================================
# Hacking output of pymatgen, because
# I'm using VASP 4.6 and not VASP 5.x
#===========================================

def get_POTCAR(poscar):
    list_postfix = ['_pv','_sv','']

    potcar = Potcar()

    list_potcar_singles = []

    for symbol in poscar.site_symbols:
        for postfix in list_postfix:
            potcar_path = PSEUDOPOTENTIALS_DIRECTORY+'%s%s/POTCAR'%(symbol,postfix)
            if os.path.isfile(potcar_path):

                single = PotcarSingle.from_file(potcar_path)

                list_potcar_singles.append(single)
                potcar.insert(-1,single)

                break

    return list_potcar_singles, potcar


def hack_poscar_file_same_element_distinct_sites(workdir,list_element_tuple):
    """
    Pymatgen is too dumb to realize that decorated elements (species) should be
    treated as different sites in VASP.

    This routine hacks the poscar file (which has already been hacked to go from 5.x to 4.6)
    to include the correct number of distinct sites.
    """
    hack_line = ''
    for (symbol, number, decoration) in list_element_tuple:
        hack_line += ' %i'%number

    with open(workdir+'POSCAR','r') as f:
        lines = f.readlines()

    lines[5] = hack_line

    with open(workdir+'POSCAR','w') as f:
        for line in lines:
            print >>f, line.strip()

    return

def hack_potcar_file(workdir,list_potcar_singles):
    """
    This function writes the potcar file again using the potcar.write_file method.
    I think VaspRun writes the potcar in the wrong order!
    """
    
    #print 're-writing the potcar'       
    #print potcar.symbols
    
    os.remove(workdir+'POTCAR')

    with open(workdir+'POTCAR','a') as f:
        for single in list_potcar_singles:
            f.write(single.data)


    return

def substitution_GGA_U_VASP(structure, workdir, Alkali='Na',nproc=16, dry_run=False, run_template=None):
    """
    This function will perform the menial tasks generating input and launching
    jobs, creating workdir along the way if it doesn't exist.

    This is more flexible than Generate_Relax_GGA_U_VASP, which will be useful
    for substituting various elements in the structure.
    """

    A = pymatgen.Element(Alkali)

    try:
        os.makedirs(workdir+'/CIF/')
    except:       
        print 'directory already exists. Continue...'

    os.chdir(workdir)

    formula  = structure.formula.replace(' ','')

    cif_filename = os.path.abspath('.')+'/CIF/'+formula+'.cif'

    structure.to(fmt='cif',filename=cif_filename)

    if not os.path.isdir(formula):
        os.mkdir(formula)
    os.chdir(formula)

    if not os.path.isdir('SGGA_U_Relaxation'):
        os.mkdir('SGGA_U_Relaxation') 
    os.chdir('SGGA_U_Relaxation') 

    file = open('run.py','w')


    print >> file, run_template.format(formula, cif_filename, nproc)
    file.close()

    if not dry_run:
        os.system('python run.py')
    os.chdir('../../../')


    return

#===========================================
# helper routines
#===========================================

def get_VASP_inputs(structure, workdir, job_name, nproc=64, kppa=500, extra_incar_dict = None):

    if os.path.exists(workdir):
        print 'WORKDIR ALREADY EXISTS. DELETE TO LAUNCH NEW JOB'
        return -1

    poscar  = Poscar(structure)

    list_potcar_singles, potcar= get_POTCAR(poscar)

    kpoints = Kpoints.automatic_density(structure, kppa=kppa)

    # Default values
    incar_dict = dict(  SYSTEM  =   structure.formula, # Name of job
                        LREAL   =   'Auto',     # Should projections be done in real space? Let VASP decide
                        ENCUT   =   520.,       # 520. eV, just like Ceder
                        IBRION  =   2,          # Controls ionic relataxion: 1-> DISS, 2 -> CG, 3-> MD
                        EDIFF   =   1E-7,       # criterion to stop SCF loop, in eV
                        EDIFFG  =  -1E-3,       # criterion to stop ionic relaxations. Negative means FORCES < |EDIFFG|
                        PREC    =   'HIGH',     # level of precision
                        AMIX    =   0.2,
                        AMIX_MAG=   0.8,
                        BMIX    =   0.001,
                        BMIX_MAG=   0.001,
                        NSW     =   150,        # Maximum number of ionic steps
                        ISMEAR  =   0,          # smearing scheme. Use 0 for insulators, as suggested by VASPWIKI
                        ISPIN   =   2,          # spin polarized 
                        NPAR    =   8,          # VASPWIKI recommends sqrt(ncore)
                        LSCALU  =   False,      # Don't use scalapack. Probably a can of worms.
                        ALGO    =   'NORMAL',   # what ionic relaxation scheme to use? 
                        LORBIT  =   11,         # 11 prints out the DOS
                        ISIF    =   3,          # Controls the computation of stress tensor. 3 computes everything
                        NSIM    =   4,          # how many bands to treat in parallel? Default is 4, probably fine.
                        SIGMA   =   0.025,      # smearing in eV
                        LMAXMIX =   4,          # Description: LMAXMIX controls up to which l-quantum number the one-center PAW charge densities are passed through the charge density mixer. MaterialsProject uses 4.
                        LCHARG  =   False,      # Write charge densities?
                        LWAVE   =   False,      # write out the wavefunctions?
                        LPLANE  =   True,       # Plane distribution of FFT coefficients. Reduces communications in FFT.
                        NELM    =   100,        # maximum number of SCF cycles.
                        NELMDL  =  -10,         # since initial orbitals may be random, fixes hamiltonian for |NELM| SCF cycles to give wf a chance to simmer down.
                        ISTART  =   0,          # begin from scratch!
                        ISYM    =   2)          # use symmetry 

    if extra_incar_dict  != None:
        incar_dict.update( extra_incar_dict  )

    incar   = Incar.from_dict(incar_dict )


    incar.write_file(workdir+'INCAR')
    poscar.write_file(workdir+'POSCAR', vasp4_compatible = True)
    kpoints.write_file(workdir+'KPOINTS')
    potcar.write_file(workdir+'POTCAR')


    potcar.sort()
    hack_potcar_file(workdir,list_potcar_singles)
    

    with open(workdir+'job.sh','w') as f:
        f.write(submit_template.format(job_name,nproc))

    with open(workdir+'clean.sh','w') as f:
        f.write(clean_template)

    return 0

def get_MaterialsProject_VASP_inputs(structure, workdir, job_name, nproc=16, supplementary_incar_dict=None):

    if os.path.exists(workdir):
        print 'WORKDIR ALREADY EXISTS. DELETE TO LAUNCH NEW JOB'
        return -1

    os.mkdir(workdir)

    input_set = MPVaspInputSet()
   
    incar  = input_set.get_incar(structure)


    # set the number of parallel processors to sqrt(nproc), 
    # as recommended in manual.
    incar.update({'NPAR':int(np.sqrt(nproc))}) 

    if supplementary_incar_dict != None:
        incar.update(supplementary_incar_dict) 

    poscar = input_set.get_poscar(structure)
    kpoints = input_set.get_kpoints(structure)
    potcar = input_set.get_potcar(structure)

    incar.write_file(workdir+'INCAR')
    poscar.write_file(workdir+'POSCAR', vasp4_compatible = True)
    kpoints.write_file(workdir+'KPOINTS')
    potcar.write_file(workdir+'POTCAR')


    with open(workdir+'job.sh','w') as f:
        f.write(submit_template.format(job_name,nproc))

    with open(workdir+'clean.sh','w') as f:
        f.write(clean_template)

    return 0

def get_ModifiedMaterialsProject_VASP_inputs(structure, workdir, job_name, nproc=16, 
                                U_strategy_instance = None, supplementary_incar_dict = None):
    """
    Inputs will be inspired by MaterialsProject, but this function is appropriate
    when we are seriously modifying the inputs such that they no longer conform to Materials Project.
    """

    if os.path.exists(workdir):
        print 'WORKDIR ALREADY EXISTS. DELETE TO LAUNCH NEW JOB'
        return -1

    os.mkdir(workdir)

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

    incar.write_file(workdir+'INCAR')
    poscar.write_file(workdir+'POSCAR', vasp4_compatible = True)
    kpoints.write_file(workdir+'KPOINTS')
    potcar.write_file(workdir+'POTCAR')


    if poscar_need_hack:
        # do we need specialized hacking of the poscar because of the U strategy?       
        new_poscar_lines = U_strategy_instance.get_new_poscar_lines()
        with open(workdir+'POSCAR','w') as f:
            for line in new_poscar_lines:
                print >>f, line.strip()

    if potcar_need_hack:
        # do we need specialized hacking of the potcar because of the U strategy?       
        new_potcar_symbols = U_strategy_instance.get_new_potcar_symbols(potcar)
        new_potcar = Potcar(new_potcar_symbols) 
        # overwrite the previous potcar
        new_potcar.write_file(workdir+'POTCAR')

    with open(workdir+'job.sh','w') as f:
        f.write(submit_template.format(job_name,nproc))

    with open(workdir+'clean.sh','w') as f:
        f.write(clean_template)

    return 0

def get_MaterialsProject_DOS_VASP_inputs(structure, previous_vasp_dir, workdir, job_name, 
                                nproc=16, kpoints_density=1000 , supplementary_incar_dict=None):

    if os.path.exists(workdir):
        print 'WORKDIR ALREADY EXISTS. DELETE TO LAUNCH NEW JOB'
        return -1
    os.mkdir(workdir)

    #================================================================================
    # The approach below leads to broken DOS, for no apparent reason
    #================================================================================
    # 
    #input_set = MPNonSCFVaspInputSet.from_previous_vasp_run(previous_vasp_dir,output_dir=workdir,
    #                                mode='Uniform', kpoints_density=kpoints_density,
    #                                user_incar_settings=supplementary_incar_dict)
    #hack_poscar_file(workdir)

    # let's use my own hack to circumvent pymatgen's weakness
    input_set = TetrahedronDosSet.from_previous_vasp_run(previous_vasp_dir,output_dir=workdir,
                                    kpoints_density=kpoints_density, user_incar_settings=supplementary_incar_dict)


    with open(workdir+'job.sh','w') as f:
        f.write(submit_template.format(job_name,nproc))

    with open(workdir+'clean.sh','w') as f:
        f.write(clean_template)

    return 0

