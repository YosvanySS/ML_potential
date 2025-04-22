import os
import re
import shutil
import subprocess

# Function to read the otput files and extract the data
def extract_data(name_file, structure_type, data_file):
    # Regular pattern to read lattice parametters
    pattern_lattice_parameter = r'lattice\s*parameter\s*\(alat\)\s*=\s*([\d\.]+)\s*a\.u\.'
    # Regular pattern to read the number of atoms
    pattern_num_atoms =  r'number of atoms/cell\s*=\s*([\d\.]+)'
    # Regular pattern to read the energy
    pattern_energy = r'total\s*energy\s*=\s*(-?\d+\.\d+)'
    # Regular pattern to read lattice vectors
    pattern_lattice_vectors = r'crystal axes: \(cart\. coord\. in units of alat\)'
    # Regular pattern to read atomic positions
    pattern_positions = 'Cartesian axes'
    # Regular pattern to read the forces
    pattern_force = r'Forces acting on atoms \(cartesian axes, Ry/au\):'
    # Regular pattern to read the stress tensor
    pattern_stress = r'total   stress'
    
    #Flags to start reading each parameter when reading the file 
    begin_extract_vectors = False
    begin_extract_positions = False
    begin_extract_forces = False
    begin_extract_stress = False
    
    #Matrix with the lattice vectors
    lattice_matrix = []
    #List with atomicspecies
    symbols = []
    #Matrix for the atomic positions
    positions = []
    #Matrix for the forces
    forces = []
    #Matrix for the stress 
    stress = []
    
    #Initializing counters for the reading
    count_positions = 0
    count_forces = 0
    count_stress = 0
    
    #Reading the outout file
    with open(name_file, 'r') as file:
        #Loop for each line in the file
        for line in file:
            #Defining coincidence for each patern with the current line
            coincidence_lattice = re.search(pattern_lattice_parameter, line)
            coincidence_num_atoms = re.search(pattern_num_atoms, line)
            coincidence_energy = re.search(pattern_energy, line)
            coincidence_lattice_vectors = re.search(pattern_lattice_vectors, line)
            coincidence_cartesian = re.search(pattern_positions, line)
            coincidence_force = re.search(pattern_force, line)
            coincidence_stress = re.search(pattern_stress, line)
            
            #If the line coincides with the pattern fot the lattice parametter
            if coincidence_lattice:
                #Reading the lattice parametter aznd converting from a. u. to anstroms
                lattice = float(coincidence_lattice.group(1))*0.52917724900001
                
            #If the line coincides with the pattern fot the number of atoms
            if coincidence_num_atoms:
                #Reading the number of atoms
                num_atoms = int(coincidence_num_atoms.group(1))
                
            #If the line coincides with the pattern fot the energy
            if coincidence_energy:
                #Reading the energy
                energy = float(coincidence_energy.group(1))
                
            #If the line coincides with the pattern fot the lattice vectors
            if coincidence_lattice_vectors:
                #Flag to start reading lattice vectors set to True
                begin_extract_vectors = True
                continue
            #If the flag to start lattice vectors extraction is true
            if begin_extract_vectors:
                #Defining the coincidence between the pattern of lines which contains lattice vectors with the current line
                coincidence_vector = re.search(r'a\(\d\)\s+=\s+\(\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)', line)
                if coincidence_vector:
                    #If there is coincidence storing the vectors and converting form alat to angstrom , then adding to a matrix
                    vector = [float(coincidence_vector.group(1))*lattice, float(coincidence_vector.group(2))*lattice, float(coincidence_vector.group(3))*lattice]
                    lattice_matrix.append(vector)
                    
            #If the line coincides with the pattern fot the atomic positions
            if coincidence_cartesian:
                #Flag to start reading atomic positions set to True
                begin_extract_positions = True
                continue
            #If the flag to start positions extraction is true
            if begin_extract_positions:
                #Defining the coincidence between the pattern of lines which contains atomic positions with the current line
                coincidence_positions = re.search(r'\s*\d+\s+(\w+)\s+tau\(\s*\d+\)\s*=\s*\(\s*([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s*\)', line)
                if coincidence_positions:
                    #If there is coincidence storing the atomic specie and  positions and converting form alat to angstrom , then adding to a matrix
                    symbol = coincidence_positions.group(1)
                    x = float(coincidence_positions.group(2))*lattice
                    y = float(coincidence_positions.group(3))*lattice
                    z = float(coincidence_positions.group(4))*lattice
                    #The counter of positions read increases
                    count_positions += 1
                    symbols.append(symbol)
                    positions.append((x, y, z))
                    #Stopping reading positions when the number of positions read get to the numlber of atoms 
                    if count_positions >= num_atoms:
                        begin_extract_positions = False
                        
            #If the line coincides with the pattern fot the forces
            if coincidence_force:
                #Flag to start reading forces set to True
                begin_extract_forces = True
                continue
            #If the flag to start force extraction is true
            if begin_extract_forces:
                #Defining the coincidence between the pattern of lines which contains the forces with the current line
                coincidence_force_components = re.search(r'atom\s+\d+\s+type\s+\d+\s+force =\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', line)
                if coincidence_force_components:
                    force_x = float(coincidence_force_components.group(1))
                    force_y = float(coincidence_force_components.group(2))
                    force_z = float(coincidence_force_components.group(3))
                    count_forces += 1
                    forces.append((force_x, force_y, force_z))
                    #Stopping reading forces when the number of forces read get to the numlber of atoms
                    if count_forces >= num_atoms:
                        begin_extract_forces = False
                        
            #If the line coincides with the pattern fot the stress
            if coincidence_stress:
                #Flag to start reading stress set to True
                begin_extract_stress = True
                continue
            #If the flag to start stress extraction is true
            if begin_extract_stress:
                #Defining the coincidence between the pattern of lines which contains the stress with the current line
                coincidence_stress_tensor = re.search(r'([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)', line)
                if coincidence_stress_tensor:
                    stress_vector = [float(coincidence_stress_tensor.group(1)), float(coincidence_stress_tensor.group(2)), float(coincidence_stress_tensor.group(3))]
                    count_stress +=1
                    stress.append(stress_vector)
                    #Stopping reading stress after 3 lines
                    if count_stress >= 3:
                        begin_extract_stress = False
    
                    
    #Writing the data
    with open(data_file, 'a') as db:
        #First line with the number of atoms
        db.write(str(num_atoms)+"\n")
        #Second line with the lattice vectors, type of configuration, energy and stress
        db.write("Lattice=\""+"{:.3f}".format(lattice_matrix[0][0])+" "+"{:.3f}".format(lattice_matrix[0][1])+" "+"{:.3f}".format(lattice_matrix[0][2])+" "+"{:.3f}".format(lattice_matrix[1][0])+" "+"{:.3f}".format(lattice_matrix[1][1])+" "+"{:.3f}".format(lattice_matrix[1][2])+" "+"{:.3f}".format(lattice_matrix[2][0])+" "+"{:.3f}".format(lattice_matrix[2][1])+" "+"{:.3f}".format(lattice_matrix[2][2])+"\" Properties=species:S:1:pos:R:3:forces:R:3 config_type="+structure_type+" energy="+str(energy)+" stress=\""+"{:.8f}".format(stress[0][0])+" "+"{:.8f}".format(stress[0][1])+" "+"{:.8f}".format(stress[0][2])+" "+"{:.8f}".format(stress[1][0])+" "+"{:.8f}".format(stress[1][1])+" "+"{:.8f}".format(stress[1][2])+" "+"{:.8f}".format(stress[2][0])+" "+"{:.8f}".format(stress[2][1])+" "+"{:.8f}".format(stress[2][2])+"\" pbc=\"T T T\"\n")
        #A line for each atom with the atomic species the position and the forces
        for i in range(num_atoms):
            db.write(''.join('{:s}'.format(symbols[i]))+"     "+"{:12.8f}".format(positions[i][0])+"     "+"{:12.8f}".format(positions[i][1])+"     "+"{:12.8f}".format(positions[i][2])+"     "+"{:12.8f}".format(forces[i][0])+"     "+"{:12.8f}".format(forces[i][1])+"     "+"{:12.8f}".format(forces[i][2])+"\n")
    return 0



#FUNCTION TO GATHER ALL THE OUTPUTS FILES IN A FOLDER TO BEGING THE DATA EXTRACTION
def gather_outputs():
    #Folder to gather all the ouputs files
    outputs_folder = 'outputs'

    # Verifying if the folder exits 
    if not os.path.exists(outputs_folder):
        # Creating the folder if it does not exits
        os.makedirs(outputs_folder)
        print(f'Directory created: {outputs_folder}')
    #If it exist removing all the files in the folder
    else:
        print(f' {outputs_folder} already exits')
        content = os.listdir(outputs_folder)
        for files in content:
            files_path = os.path.join(outputs_folder, files)
            if os.path.isfile(files_path):
                os.remove(files_path)
    main_folder_path = os.getcwd()
    output_path = os.path.join(main_folder_path, outputs_folder)


    # Loop for all the folders in the main folder to find the outputs
    for folder in os.listdir(main_folder_path):
        folder_path = os.path.join(main_folder_path, folder)
        if os.path.isdir(folder_path):
            # Path to the folder already_done in the current directory
            path_already_done = os.path.join(folder_path, "already_done")
            if os.path.exists(path_already_done):
                # Loop for all the files in already_done directory
                for file in os.listdir(path_already_done):
                    #If it is an .out file the cppying to outputs directory
                    if file.endswith(".out"):
                    # Copiar el archivo a la carpeta actual
                        shutil.copy(os.path.join(path_already_done, file), output_path)
    return outputs_folder




#FUNCTION TO CREATE A LIST WITH ALL THE OUTPUT FILES NAMES
def create_out_list(directory):   
    type_calc = []
    
    # Walk through the directory for crating a txt file with the list of all outputs
    for path, subdirs, files in os.walk(directory):
        for filename in files:
            full_name_splitted = filename.split('_')
            val_type_calc = '_'.join(full_name_splitted[:-1])
            type_calc.append(val_type_calc)
    
    return type_calc, files


##########################     MAIN CODE     ################################

#Name of the file to write the data
data_file = 'train.xyz'

# Deleting the file if it exists
# if os.path.exists(data_file):
#     os.remove(data_file)

# Gathers all the ouputs files in a directory
# directory = gather_outputs()

directory = './outputs'

os.remove('./outputs/train.xyz')
os.remove('./outputs/list_out.txt')

# Creating a list with all the ouputs names
type_calc, out_files = create_out_list(directory)

os.chdir('./outputs')
# Iteration for all the outputs files to extract the data

list_out = open("list_out.txt", "a")

for i, file in enumerate(out_files):
    with open("list_out.txt", "a") as filetxt:
        line = str(i) + ' ' + file + '\n'
        filetxt.write(line)
    extract_data(file, type_calc[i], data_file)








