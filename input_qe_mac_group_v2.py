#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import os
import subprocess as sbp
import time
import shutil
import linecache as lnch
import numpy as np
import math as mt
import psutil

'-----------------------------------'
'------------- CLASSES -------------'
'-----------------------------------'
class QEWriter:
    '-----------------------------------------------------------------------'
    # Initializing the instances
    def __init__(self, qe_file_name, nat, vect_matr, atm_names, pos_matr_ang, machine, pseudo_cluster, pseudo_local, forc_conv_thr, nstep, nspin, ecutwfc, ecutrho, conv_thr, mixing_beta, dk, units):
        self.qe_file_name   = qe_file_name
        self.nat            = nat
        self.vect_matr      = vect_matr
        self.atm_names      = atm_names
        self.pos_matr_ang   = pos_matr_ang
        self.machine        = machine
        self.pseudo_cluster = pseudo_cluster
        self.pseudo_local   = pseudo_local
        self.forc_conv_thr  = forc_conv_thr
        self.nstep          = nstep
        self.nspin          = nspin
        self.ecutwfc        = ecutwfc
        self.ecutrho        = ecutrho
        self.conv_thr       = conv_thr
        self.mixing_beta    = mixing_beta
        self.dk             = dk
        self.units          = units
    
    
    
    '-----------------------------------------------------------------------'
    # Method to convert atomic coordinates from Angstrom to crystal units
    def ang2crystal(self):
        pos_matr_crys = np.zeros((3, self.nat))                         # Declaring the matrix with atomic positions in crystal units
        vect_matr_inv = np.linalg.inv(self.vect_matr)                   # Inverse matrix of the vector matrix
        pos_matr_crys = np.dot(vect_matr_inv, self.pos_matr_ang)        # Converting Angstrom positions to crystal units
        
        return pos_matr_crys
    
    
    
    '-----------------------------------------------------------------------'
    # Method for determining k-points grid fixing a delta_k between the points
    def kpoints(self):
        # vector with k-points
        k_points = np.zeros(3)
        # Determining the values of k-points in X, Y and Z directions depending of the value of delta_k
        for i in range(3):
            norm_vect = np.sqrt(np.sum(np.power(self.vect_matr[:, i], 2)))      # Norm of the lattice vectors in real space: sqrt(a_x^2 + a_y^2 + a_z^2), same for b and c vectors
            norm_vect_kspace = 2*mt.pi / norm_vect                              # The norm of the lattice vectors in the reciprocal space
            k_points[i] = mt.ceil(norm_vect_kspace / self.dk)                   # This give us the number of k-points depending of the cosidered delta_k. The value is rounded to the next integer.
        
        return k_points
    
    
    
    '-----------------------------------------------------------------------'
    # Method to determine the number of electonic bands nbnd
    def number_bnd(self):
        # Variable to count the number of valance and core electrons
        _nelec_ = 0
        # Determining the total number of valance and core electrons
        for atm_type in self.atm_names:
            if atm_type == 'W':
                _nelec_ += 14                   # 14 electrons considered in W atoms
            elif atm_type == 'Cu':
                _nelec_ += 11                   # 11 electrons considered in Cu atoms
            elif atm_type == 'H':
                _nelec_ += 1                    # 1 electron considered in Cu atoms
        # Determining the number of elctronic bands, which is equal to the number of electrons/2 + a 20% of nbnd for metals regarding pw.x webpage
        nbnd = mt.ceil(0.6 * _nelec_)
        
        return nbnd



    '-----------------------------------------------------------------------'
    # Method to define how many type of atoms are in the structure and save the pseudopential for QE files
    def number_typ(self):
        # Set the similar atoms names in unique_symbols list
        unique_symbols = set(self.atm_names)
        ntyp = len(unique_symbols)
        pseudos = ['W   183.84 W.pbe-nsp-van.UPF', 'Cu  63.546 Cu.pbe-dn-rrkjus_psl.1.0.0.UPF']
        pseudo_qe = []
        
        if ntyp == 1:
            if self.atm_names[0] == 'W':
                pseudo_qe.append(pseudos[0])
            elif self.atm_names[0] == 'Cu':
                pseudo_qe.append(pseudos[1])
        else:
            pseudo_qe = pseudos
        
        return ntyp, pseudo_qe
    
    

    '-----------------------------------------------------------------------'
    # Method to write the input file of QE
    def writing_input(self):
        # Saving outputs of methods
        ntyp, pseudo_qe = self.number_typ()
        nbnd = self.number_bnd()
        k_points = self.kpoints()
        k_points = k_points.astype(int)             # Transform k_points in integer type
        # Saving the paths to the pseudoptentials for the cluster and the local machine
        if self.machine == 'cluster':
            folder_path = './cluster'
            pseudo_dir = self.pseudo_cluster
        else:
            folder_path = './local'
            pseudo_dir = self.pseudo_local
        # Joinging the path of the folder to the name of the input file to create the file
        file_path = os.path.join(folder_path, self.qe_file_name)
        
        # Open the file to write the data in it
        with open(file_path+'.in', 'w') as file_qe:
            # Writing &control parameters
            file_qe.write("&control\n")
            file_qe.write("    calculation   = 'scf'\n")
            file_qe.write("    prefix        = '"+str(self.qe_file_name)+"'\n")
            file_qe.write("    pseudo_dir    = '"+pseudo_dir+"'\n")
            file_qe.write("    verbosity     = 'high'\n")
            file_qe.write("    restart_mode  = 'from_scratch'\n")
            file_qe.write("    forc_conv_thr = "+str("{:.2e}".format(self.forc_conv_thr))+",\n")
            file_qe.write("    wf_collect    = .true.\n")
            file_qe.write("    tstress       = .true.\n")
            file_qe.write("    tprnfor       = .true.\n")
            file_qe.write("    nstep         = "+str(self.nstep)+"\n")
            file_qe.write("/\n")
            
            # Writing &system parameters
            file_qe.write("&system\n")
            file_qe.write("    ibrav = 0,\n")
            file_qe.write("    nat = "+str(self.nat)+",\n")
            file_qe.write("    ntyp = "+str(ntyp)+",\n")
            file_qe.write("    nbnd = "+str(nbnd)+",\n")
            file_qe.write("    occupations = 'smearing',\n")
            file_qe.write("    degauss = 2d-2,\n")
            file_qe.write("    smearing = 'mv',\n")
            file_qe.write("    nspin = "+str(self.nspin)+",\n")
            file_qe.write("    ecutwfc = "+str(self.ecutwfc)+"\n")
            file_qe.write("    ecutrho = "+str(self.ecutrho)+"\n")
            file_qe.write("/\n")
            
            # Writing &electrons parameters
            file_qe.write("&electrons\n")
            file_qe.write("    diagonalization = 'cg',\n")
            file_qe.write("    conv_thr = "+str("{:.2e}".format(self.conv_thr))+"\n")
            file_qe.write("    mixing_beta = "+str(self.mixing_beta)+"\n")
            file_qe.write("/\n")
            
            # Writing &ions parameters
            file_qe.write("&ions\n")
            file_qe.write("    ion_dynamics = 'bfgs'\n")
            file_qe.write("/\n")
            
            # Writing ATOMIC SPECIES parameters
            file_qe.write("\n")
            file_qe.write("ATOMIC_SPECIES\n")
            if ntyp == 1:
                file_qe.write(pseudo_qe[0]+"\n")
            elif ntyp == 2:
                for i in range(ntyp):
                    file_qe.write(pseudo_qe[i]+"\n")
            file_qe.write("\n")
            
            # Writing K_POINTS AUTOMATIC parameters
            file_qe.write("K_POINTS AUTOMATIC\n")
            file_qe.write(""+str(k_points[0])+" "+str(k_points[1])+" "+str(k_points[2])+" 1 1 1\n")
            file_qe.write("\n")
            
            # Writing CELL_PARAMETERS parameters
            file_qe.write("CELL_PARAMETERS (angstrom)\n")
            for i in range(3):
                file_qe.write("{0:14.10f} {1:14.10f} {2:14.10f}\n"
                            .format(self.vect_matr[0, i], self.vect_matr[1, i], self.vect_matr[2, i]))
            file_qe.write("\n")
            
            # Writing ATOMIC_POSITIONS parameters
            file_qe.write("ATOMIC_POSITIONS ("+self.units+")\n")
            if self.units == 'crystal':
                cryst_unit = self.ang2crystal()
                for i in range(self.nat):
                    file_qe.write("{0:2} {1:14.10f} {2:14.10f} {3:14.10f}\n"
                                .format(self.atm_names[i], cryst_unit[0, i], cryst_unit[1, i], cryst_unit[2, i]))
            elif self.units == 'angstrom':
                for i in range(self.nat):
                    file_qe.write("{0:2} {1:14.10f} {2:14.10f} {3:14.10f}\n"
                                .format(self.atm_names[i], self.pos_matr_ang[0, i], self.pos_matr_ang[1, i], self.pos_matr_ang[2, i]))
        
        file_qe.close()
        
        return 0





'==========================================================================='
class FileReader:
    '-----------------------------------------------------------------------'
    # Initializing the instances
    def __init__(self, file_name, pseudo_cluster, pseudo_local, forc_conv_thr, nstep, nspin, ecutwfc, ecutrho, conv_thr, mixing_beta, dk, units):
        self.file_name      = file_name
        self.pseudo_cluster = pseudo_cluster
        self.pseudo_local   = pseudo_local
        self.forc_conv_thr  = forc_conv_thr
        self.nstep          = nstep
        self.nspin          = nspin
        self.ecutwfc        = ecutwfc
        self.ecutrho        = ecutrho
        self.conv_thr       = conv_thr
        self.mixing_beta    = mixing_beta
        self.dk             = dk
        self.units          = units



    '-----------------------------------------------------------------------'
    # Method to read the file with the db
    def read_file(self):
        
        # lists to save the cell primitive vectors from the data base (db)
        a, b, c = [], [], []
        
        # Declaring the matrix of the cell primitive vectors: a, b, c. Those values are read form the db.
        #[a_x b_x c_x
        # a_y b_y c_y
        # a_z b_z c_z]
        vect_matr = np.zeros((3, 3))
        type_calc = []
        
        # Open the file to count the number of structures: cont_str
        with open(self.file_name, 'r') as file:
            data = file.read()                  # Data inside the file
            find = 'Lattice'                    # String to find into the file
            cont_str = data.count(find)         # Counting the number of times that find appears, which is equals to the number of structures.
        file.close()
        
        with open(self.file_name, 'r') as file:
            for i in range(cont_str):
                line1 = file.readline()         # Reading the number of atoms in each structure: nat
                nat = int(line1)
                
                # Reading line 2 of each strucuture to save the values of the cell primitive vectors: a, b, c.
                line2 = file.readline()
                line_split = re.split('[=" ]', line2)       # Splitting line2 using [=" ] separators
                
                # Creating the cell primitive vector matrix: vect_matr
                for n in range(3):
                    for m in range(3):
                        vect_matr[m, n] = float(line_split[3*n+m+2])
                
               # Determine the position of the calculation type in line 2
                position_type_calc=line_split.index("config_type")+1    
                type_calc.append(line_split[position_type_calc])
                
                pos_matr_ang = np.zeros((3, nat))   # Matrix with atomic positions in Angstorm.
                atm_names = []                      # Row vector with the atomic species read from the file.
                
                # Reading the atomic coordinates
                for j in range(nat):
                    line = file.readline()
                    line_split = line.split()
                    
                    atm_names.append(line_split[0])     # Saving the atomic species in atm_names
                    
                    for k in range(3):
                        pos_matr_ang[k, j] = float(line_split[k+1])     # Saving atomic positions as column vectors
                        
                qe_name = type_calc[i]+'_'+str(i+1)
                
                writer = QEWriter(qe_name, nat, vect_matr, atm_names, pos_matr_ang, 'cluster', self.pseudo_cluster, self.pseudo_local, self.forc_conv_thr, self.nstep, self.nspin, self.ecutwfc, self.ecutrho, self.conv_thr, self.mixing_beta, self.dk, self.units)
                writer.writing_input()
                writer = QEWriter(qe_name, nat, vect_matr, atm_names, pos_matr_ang, 'local', self.pseudo_cluster, self.pseudo_local, self.forc_conv_thr, self.nstep, self.nspin, self.ecutwfc, self.ecutrho, self.conv_thr, self.mixing_beta, self.dk, self.units)
                writer.writing_input()
        
        file.close()
        
        return cont_str, type_calc





'==========================================================================='
class RunCreator:
    '-----------------------------------------------------------------------'
    # Initializing the instances
    def __init__(self, cont_str, type_calc, time_exe, folder1_path, folder2_path, path_qe):
        self.cont_str     = cont_str
        self.type_calc    = type_calc
        self.time_exe     = time_exe
        self.folder1_path = folder1_path
        self.folder2_path = folder2_path
        self.path_qe      = path_qe
    
    
    
    '-----------------------------------------------------------------------'
    # Method to create run files already parallelized
    def find_divisors(self, val):
        divisors = []
        
        # Searching from 1 to the square root of the number
        for i in range(1, int(mt.sqrt(val)) + 1):
            if val % i == 0:
                divisors.append(i)
                # If 'i' is not the square root of 'number', add the corresponding divisor
                if i != val // i:
                    divisors.append(val // i)
        
        divisors.sort()
        
        return divisors

    
    
    '-----------------------------------------------------------------------'
    # Method to execute QE to determine parallelization parameters
    def pwx_qe(self):
        os.chdir(self.folder2_path)       # Changing the repository to execute pw.x code in the local repository        
        
        for i in range(self.cont_str):
            input_i = self.type_calc[i]+'_'+str(i+1) + '.in'
            output_i = self.type_calc[i]+'_'+str(i+1) + '.out'
            command = ""+self.path_qe+"pw.x <"+input_i+" >"+output_i+""
            
            process = sbp.Popen(command, shell=True)
            
            # Sleeping the process during 2s of execution
            time.sleep(2)
            
            # Checking if the process is still running
            if process.poll() is None:
                process.kill()
                print('Already generated:', output_i)
        
        os.chdir('../')
            
        return 0
    
    
    
    '-----------------------------------------------------------------------'
    # Method to write the run files for IDRIS already parallelized in the create_run method
    def run_writer_idris(self, name_job, npools, ntasks, nodes):
        run_name = 'run-'+name_job
        
        with open(run_name, 'w') as file_run:
            # Writing lines of run file
            file_run.write("#!/bin/bash\n")
            file_run.write("\n")
            
            file_run.write("#SBATCH --account=mzn@cpu\n")
            file_run.write("#SBATCH --job-name="+name_job+"\n")
            file_run.write("#SBATCH --output="+name_job+".out           # nom du fichier de sortie\n")
            file_run.write("#SBATCH --error="+name_job+"_%j.err          # nom du fichier d'erreur (ici en commun avec la sortie)\n")
            file_run.write("\n")
            
            file_run.write("#SBATCH --ntasks="+str(ntasks)+"\n")
            file_run.write("#SBATCH --ntasks-per-node=40\n")
            file_run.write("#SBATCH --hint=nomultithread\n")
            file_run.write("#SBATCH --time="+str(self.time_exe)+"\n")
            file_run.write("#SBATCH --qos=qos_cpu-t3\n")
            file_run.write("\n")
            
            file_run.write("cd ${SLURM_SUBMIT_DIR}\n")
            file_run.write("\n")
            
            file_run.write("module purge\n")
            file_run.write("module load intel-compilers\n")
            file_run.write("module load intel-mpi\n")
            file_run.write("module load quantum-espresso/6.4-mpi\n")
            file_run.write("\n")
            
            file_run.write("set -x\n")
            file_run.write("\n")
            
            file_run.write("DIR=$(pwd)\n")
            file_run.write("\n")
            
            file_run.write("mkdir $WORK/ML/\n")
            file_run.write("cp "+name_job+".in $WORK/ML/\n")
            file_run.write("cd $WORK/ML/\n")
            file_run.write("\n")
            
            file_run.write("srun -n "+str(ntasks)+" pw.x -npools "+str(npools)+" -input "+name_job+".in > "+name_job+".out\n")
            file_run.write("\n")
            file_run.write("cp "+name_job+".out $DIR\n")
            
        file_run.close()
        
        #os.chdir('../')
        
        return 0
    
    # Method to write the run files for ADASTRA already parallelized in the create_run method
    def run_writer_adastra(self, name_job, npools, ntasks, nodes):
        run_name = 'run-'+name_job
        
        with open(run_name, 'w') as file_run:
            # Writing lines of run file
            file_run.write("#!/bin/bash\n")
            file_run.write("\n")
            file_run.write("#SBATCH --account=c1606612\n")
            file_run.write("#SBATCH --job-name="+name_job+"\n")
            file_run.write("#SBATCH --constraint=GENOA\n")
            file_run.write("#SBATCH --exclusive\n")
            file_run.write("#SBATCH --ntasks="+str(ntasks)+"\n")
            file_run.write("#SBATCH --ntasks-per-node=192\n")
            file_run.write("#SBATCH --cpus-per-task=1\n")
            file_run.write("#SBATCH --threads-per-core=1\n")
            file_run.write("#SBATCH --error="+name_job+"_%j.err          # nom du fichier d'erreur (ici en commun avec la sortie)\n")
            file_run.write("#SBATCH --output="+name_job+".out           # nom du fichier de sortie\n")
            file_run.write("#SBATCH --time="+str(self.time_exe)+"\n")
            file_run.write("\n")
            
            
            file_run.write("module purge\n")
            file_run.write("\n")
            
            file_run.write("module load develop\n")
            file_run.write("module load CCE-CPU-3.1.0\n")
            file_run.write("module load PrgEnv-intel/8.5.0\n")
            file_run.write("\n")
            
            file_run.write("set -x\n")
            file_run.write("\n")
            
            file_run.write("DIR=$(pwd)\n")
            file_run.write("TARGET="+name_job+"\n")
            file_run.write("BINARY=/lus/home/CT8/c1606612/yferro/QE/qe-6.8/bin/pw.x\n")
            file_run.write("\n")
            
            file_run.write("rm -rf $WORKDIR/$TARGET\n")
            file_run.write("mkdir $WORKDIR/$TARGET\n")
            file_run.write("cp $DIR/${TARGET}.in $WORKDIR/$TARGET\n")
            file_run.write("cd $WORKDIR/$TARGET\n")
            file_run.write("\n")
            
            file_run.write("export OMP_PROC_BIND=CLOSE\n")
            file_run.write("export OMP_PLACES=THREADS\n")
            file_run.write("\n")
            
            file_run.write("srun -- $BINARY -npools "+str(npools)+" -input ${TARGET}.in > ${TARGET}.out\n")
            file_run.write("\n")
            file_run.write("cp ${TARGET}.out $DIR\n")
            file_run.write("done_flag=$(grep 'JOB DONE' ${TARGET}.out | wc -l)\n")
            file_run.write("conv_flag=$(grep 'convergence NOT achieved after' ${TARGET}.out | wc -l)\n")
            file_run.write("stop_flag=$(grep 'Program stopped by user request' ${TARGET}.out | wc -l)\n")
            file_run.write("cd $WORKDIR\n")
            file_run.write("\t if [[ ! $TARGET =~ ^phonon ]]; then\n")
            file_run.write("\t \t if [ \"$done_flag\" -eq 1 ] && [ \"$conv_flag\" -eq 0 ] && [ \"$stop_flag\" -eq 0 ]; then\n")
            file_run.write("\t \t \t rm -rf $WORKDIR/$TARGET\n")
            file_run.write("\t \t \t echo \"Folder '$TARGET' erased.\"\n")
            file_run.write("\t \t else\n")
            file_run.write("\t \t \t mv $WORKDIR/$TARGET $SCRATCHDIR\n")
            file_run.write("\t \t \t echo \" '$TARGET' not converged\"\n")
            file_run.write("\t \t fi\n")
            file_run.write("\t else\n")
            file_run.write("\t \t mv $TARGET $SCRATCHDIR\n")
            file_run.write("\t \t echo \"Folder '$TARGET' not erased 'phonons'. Moved to $SCRATCHDIR\"\n")
            file_run.write("\t fi\n")
            
        file_run.close()
        return 0
    
    
    
    '-----------------------------------------------------------------------'
    # Method to write the run files already parallelized in the create_run method
    def run_writer_general(self, npools, ntasks, nodes):
        run_name = 'run-ntasks_'+str(ntasks)
        
        with open(run_name, 'w') as file_run:
            # Writing lines of run file
            file_run.write("#!/bin/bash\n")
            file_run.write("\n")
            file_run.write("#SBATCH --account=c1606612\n")
            file_run.write("#SBATCH --job-name=nt_"+str(ntasks)+"\n")
            file_run.write("#SBATCH --constraint=GENOA\n")
            file_run.write("#SBATCH --exclusive\n")
            file_run.write("#SBATCH --ntasks="+str(ntasks)+"\n")
            file_run.write("#SBATCH --ntasks-per-node=192\n")
            file_run.write("#SBATCH --cpus-per-task=1\n")
            file_run.write("#SBATCH --threads-per-core=1\n")
            file_run.write("#SBATCH --error=nt_"+str(ntasks)+"_%j.err          # nom du fichier d'erreur (ici en commun avec la sortie)\n")
            file_run.write("#SBATCH --output=nt_"+str(ntasks)+".out           # nom du fichier de sortie\n")
            file_run.write("#SBATCH --time="+str(self.time_exe)+"\n")
            file_run.write("\n")
            
            
            file_run.write("module purge\n")
            file_run.write("\n")
            
            file_run.write("module load develop\n")
            file_run.write("module load CCE-CPU-3.1.0\n")
            file_run.write("module load PrgEnv-intel/8.5.0\n")
            file_run.write("\n")
            
            file_run.write("set -x\n")
            file_run.write("\n")
            
            file_run.write("DIR=$(pwd)\n")
            file_run.write("BINARY=/lus/home/CT8/c1606612/yferro/QE/qe-6.8/bin/pw.x\n")
            file_run.write("\n")
            
            file_run.write("rm already_done.dat\n")
            file_run.write("mkdir -p already_done\n")
            file_run.write("\n")
            
            file_run.write("input_files=\"input_files.dat\"\n")
            file_run.write("name_file=\"name_file.dat\"\n")
            file_run.write("already_done=\"already_done.dat\"\n")
            file_run.write("pending=\"pending.dat\"\n")
            file_run.write("\n")
            
            file_run.write("ls *.in > \"$input_files\"\n")
            file_run.write("sed 's/...$//' \"$input_files\" > \"$name_file\"\n")
            file_run.write("\n")
            
            file_run.write("touch \"$already_done\"\n")
            file_run.write("touch \"$pending\"\n")
            file_run.write("\n")
            
            file_run.write("i=0\n")
            file_run.write("while IFS= read -r name_i; do\n")
            file_run.write("\t i=$((i+1))\n")
            file_run.write("\t name_list[$i]=\"$name_i\"\n")
            file_run.write("done < \"$name_file\"\n")
            file_run.write("\n")
            
            file_run.write("rm \"$name_file\"\n")
            file_run.write("rm \"$input_files\"\n")
            file_run.write("\n")
            
            file_run.write("for j in $(seq 1 $i); do\n")
            file_run.write("\t rm \"$pending\"\n")
            file_run.write("\n")
            
            file_run.write("\t for ((k = j; k <= i; k++)); do\n")
            file_run.write("\t \t echo \"${name_list[k]}\" >> $pending\n")
            file_run.write("\t done\n")
            file_run.write("\n")
            
            file_run.write("\t TARGET=\"${name_list[j]}\"\n")
            file_run.write("\t rm -rf $WORKDIR/$TARGET\n")
            file_run.write("\t mkdir $WORKDIR/$TARGET\n")
            file_run.write("\t cp $DIR/${TARGET}.in $WORKDIR/$TARGET\n")
            file_run.write("\t cd $WORKDIR/$TARGET\n")
            file_run.write("\n")
            
            file_run.write("\t srun -- $BINARY -npools "+str(npools)+" -input ${TARGET}.in > ${TARGET}.out\n")
            file_run.write("\n")
            
            file_run.write("\t cp ${TARGET}.out $DIR\n")
            file_run.write("\n")
            
            file_run.write("\t done_flag=$(grep 'JOB DONE' ${TARGET}.out | wc -l)\n")
            file_run.write("\t conv_flag=$(grep 'convergence NOT achieved after' ${TARGET}.out | wc -l)\n")
            # file_run.write("stop_flag=$(grep 'Program stopped by user request' ${TARGET}.out | wc -l)\n")
            file_run.write("\n")
            
            file_run.write("\t cd $WORKDIR\n")
            file_run.write("\n")
            
            # file_run.write("\t if [[ ! $TARGET =~ ^phonon ]]; then\n")
            file_run.write("\t if [ \"$done_flag\" -eq 1 ] && [ \"$conv_flag\" -eq 0 ] && [ \"$stop_flag\" -eq 0 ]; then\n")
            file_run.write("\t \t rm -rf $WORKDIR/$TARGET\n")
            file_run.write("\t \t echo \"Folder '$TARGET' erased.\"\n")
            file_run.write("\t else\n")
            file_run.write("\t \t mv $WORKDIR/$TARGET $SCRATCHDIR\n")
            file_run.write("\t \t echo \" '$TARGET' not converged\"\n")
            file_run.write("\t fi\n")
            # file_run.write("\t else\n")
            # file_run.write("\t \t mv $TARGET $SCRATCHDIR\n")
            # file_run.write("\t \t echo \"Folder '$TARGET' not erased 'phonons'. Moved to $SCRATCHDIR\"\n")
            # file_run.write("\t fi\n")
            file_run.write("\n")
            
            file_run.write("\t cd $DIR\n")
            file_run.write("\n")
            
            file_run.write("\t mv \"${name_list[j]}\".* already_done/\n")
            file_run.write("\t cd already_done/\n")
            file_run.write("\t ls *.in > temp.dat\n")
            file_run.write("\t sed 's/...$//' temp.dat > \"$already_done\"\n")
            file_run.write("\t rm temp.dat\n")
            file_run.write("\t mv \"$already_done\" ../\n")
            file_run.write("\t cd ..\n")
            file_run.write("done\n")
            file_run.write("rm \"$pending\"\n")
            file_run.write("echo 'No jobs pending' > $pending\n")
            
        file_run.close()
        
        #os.chdir('../')
        
        return 0



    '-----------------------------------------------------------------------'
    # Method to create run files already parallelized
    def create_run(self):
        
        # # Folders to subgroup the files according to the same ntasks
        group_files = []
        
        # Pattertns to search in output files for the parallelization
        pattern1 = 'number of k points='
        pattern2 = 'Dense  grid:'
        ntasks_per_node = 40
            
        # Rerading every output file
        for i in range(self.cont_str):
            output_i = self.type_calc[i] + '_' + str(i+1) + '.out'
            run_i    = 'run-' + str(i+1)
            with open(output_i, 'r') as file:
                # Reading every line in the output file to find the pattern1 and pattern2
                for line in file:
                    # If pattern1 is found in line
                    if re.search(pattern1, line):
                        line_split = line.split()
                        npools = int(line_split[4])             # Number of polls = number of k-points
                    if re.search(pattern2, line):
                        line_split = line.split()
                        FFT_z = line_split[9]                   # Number of Fast Fourier Transforms in z direction
                        FFT_z = int(FFT_z.rstrip(')'))          # Removing the bracket which FFT_z is saved from the output file
                        
                # ngroups is the number of groups that we parallelize the FFT in the z direction. This number has to be a divisor of FFT_z
                FFT_z_div = self.find_divisors(FFT_z)
                
                if npools >= 300:
                    npools = 300
                
                # Backward loop for detecting the highest divisor that can be used
                for j in range(len(FFT_z_div)-1, -1, -1):
                    if FFT_z_div[j]*npools < 3000:
                        ngroups = FFT_z_div[j]
                        break
                    else:
                        continue
                ntasks = ngroups * npools
                ncpus = ntasks
                nodes = mt.ceil(ncpus / ntasks_per_node)
                
                sub_folder_name = 'ntasks_'+str(ntasks)
                moving_input = self.type_calc[i]+'_'+str(i+1)+'.in'
                
                # If the ntasks is not in group_files create a directory called ntasks
                if ntasks not in group_files:
                    group_files.append(ntasks)
                    os.mkdir(sub_folder_name)
                    os.chdir(sub_folder_name)
                    self.run_writer_general(npools, ntasks, nodes)
                else:
                    os.chdir(sub_folder_name)
                name_job = self.type_calc[i]+'_'+str(i+1)    
                
                self.run_writer_adastra(name_job, npools, ntasks, nodes)
                sbp.run('cp ../../cluster/'+moving_input+' ./', shell=True)
                
                os.chdir('../')

        file.close()
        
        moving_run = '../'+folder1_path
        sbp.run('mv ntasks_* '+moving_run, shell=True)
        os.chdir('../')
                
        return 0





'==========================================================================='
class BashCreator:
    '-----------------------------------------------------------------------'
    # Initializing the instances
    def __init__(self, db_file, bash_file, path_qe, path_parallel, path_timeout, ncpus, type_calc):
        self.db_file   = db_file
        self.bash_file = bash_file
        self.path_qe   = path_qe
        self.path_parallel=path_parallel
        self.path_timeout=path_timeout
        self.ncpus     = ncpus
        self.type_calc = type_calc
    
    
    
    '-----------------------------------------------------------------------'
    # Method to create bash file to compile QE from local machine
    def temp_bash(self):
        bash_name = self.bash_file
        
        with open("name.dat", 'w') as file_temp:
            for i, name in enumerate(self.type_calc):
                file_temp.write(name+'_'+str(i+1)+'\n')
        
        with open(bash_name, 'w') as file_bash:
            
            # Writing lines of bash file
            file_bash.write("#!/bin/bash\n")
            file_bash.write("\n")
            
            #Determine the number of configurations
            file_bash.write("nconf=$(grep -o 'Lattice=' ../"+self.db_file+" | wc -l)\n")
            file_bash.write("\n")
            
            # Creating the 'run_and_echo' function in bash to run QE during 2 seconds and print in console the creation of output files.
            file_bash.write("run_and_echo() {\n")
            file_bash.write(""+self.path_timeout+"timeout 2s "+self.path_qe+"pw.x < $1.in > $1.out; echo \"File $1.out generated\"\n")
            file_bash.write("}\n")
            file_bash.write("\n")
            
            # Exporting 'run_and_echo'
            file_bash.write("export -f run_and_echo\n")
            file_bash.write("\n")
            
            file_bash.write("i=0\n")
            file_bash.write("\n")
            
            #Reading the file names from name.dat
            file_bash.write("declare -a name_file\n")
            file_bash.write("while IFS= read -r LINE;\n")
            file_bash.write("do\n")
            file_bash.write("i=$((i+1))\n")
            file_bash.write("name_file[$i]=\"$LINE\"\n")
            file_bash.write("done < name.dat\n")
            
            #Calling the function 'run and echo' for each input file
            #file_bash.write("for i in $(seq 1 $nconf)\n")
            #file_bash.write("do\n")
            file_bash.write("printf \"%s\n\" \"${name_file[@]}\" | "+self.path_parallel+"parallel --will-cite --keep-order -j"+str(self.ncpus)+" run_and_echo 2>/dev/null\n")
            #file_bash.write("done\n")
            
        file_bash.close()
        
        return 0





'==========================================================================='
class Sorting:
    '-----------------------------------------------------------------------'
    # Initializing the instances
    def __init__(self, folder_path, cont_str):
        self.folder_path = folder_path
        self.cont_str    = cont_str
    
    
    
    # '-----------------------------------------------------------------------'
    # # Method to create bash file to compile QE from local machine
    # def classing(self):
        
    






'==========================================================================='
'-------------------------------- MAIN CODE --------------------------------'
'==========================================================================='
if __name__ == "__main__":
    # Folders for saving files
    folder1_path = './cluster'
    folder2_path = './local'
    
    paths = [folder1_path, folder2_path]
    
    for folder in paths:
        if not os.path.exists(folder):
            # If the folder doesn't exist
            os.makedirs(folder)
            print('Folders were created correctly!')
        else:
            print('Folders already exits!')
    
    # File with db
    db_file = 'train.xyz'
    
    # Bash file to execute QE
    bash_file = 'qe_exe.sh'
    
    # Number of cpus to RUN QE IN LOCAL MACHINE
    ncpus = 8
    
    
    #-----------------------------------------------------------------------#
    #------------------------- QE input variables --------------------------#
    #-----------------------------------------------------------------------#
    # paths to the psudopotential folder: pseudo_cluster and pseudo_local for the cluster and local machine respectively
    #pseudo_cluster = '/linkhome/rech/genpii01/uex56lt/Pseudo' #IDRIS
    pseudo_cluster = '/lus/home/CT8/c1606612/yferro/Pseudo'    #ADASTRA
    pseudo_local   = '/Users/yosvany/Documents/Yosvany/PhD/HPC/Pseudo'
    
    # Path to the executable of QE in the local machine
    #path_qe = '/opt/Softs/QE/qe-6.8/bin/'
    path_qe='/opt/Softs/QE/qe-6.8/bin/'
    path_parallel='/opt/local/bin/'
    path_timeout='/opt/local/bin/'
    
    # &control parameters
    forc_conv_thr = 1.0e-4
    nstep = 150
    
    # &system parameters
    nspin = 1
    ecutwfc = 50.0
    ecutrho = 400.0
    
    # &electrons parameters
    conv_thr = 1.0e-6
    mixing_beta = 0.2
    
    # K_POINTS parameters
    dk_cu = 0.072
    
    # Atomic units system
    units = 'crystal'
    
    # Time execution of QE in cluster. Variable for the run files. Format hh:mm:ss
    time_exe = '20:00:00'
    
    #-----------------------------------------------------------------------#
    #--------------------------- Reading database --------------------------#
    #-----------------------------------------------------------------------#
    reader = FileReader(db_file, pseudo_cluster, pseudo_local, forc_conv_thr, nstep, nspin, ecutwfc, ecutrho, conv_thr, mixing_beta, dk_cu, units)
    cont_str, type_calc = reader.read_file()
    
    print('----------------------------------')
    print('Number of strucutres:', cont_str)
    print('----------------------------------')
    

    #-----------------------------------------------------------------------#
    #------------------------- Creating bash script ------------------------#
    #-----------------------------------------------------------------------#    
    # Changing the repository to create and execute the bash script in the local repository
    os.chdir(folder2_path) 

    bash = BashCreator(db_file, bash_file, path_qe, path_parallel, path_timeout, ncpus, type_calc)
    bash.temp_bash()
     
    # Executing command in console
    command = "bash "+bash_file+""
    
    process = sbp.run(command, shell=True)
    if process.returncode == 0:
        print("The outputs were correctly generated")
    else:
        print("Error generating the outputs")
    
    
    #-----------------------------------------------------------------------#
    #-------------------------- Creating run files -------------------------#
    #-----------------------------------------------------------------------# 
    runner = RunCreator(cont_str, type_calc, time_exe, folder1_path, folder2_path, path_qe)
    runner.create_run()
    print('----------------------------------------------------')
    print('Run files created correctly')
    print('----------------------------------------------------')
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
