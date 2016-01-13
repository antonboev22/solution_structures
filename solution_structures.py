#! /usr/bin/env python3
'''
\section{solution\_structures.py}

\\hyperlink{content}{Content}

Program is intended for making of the solid solution from some initial lattice

Requirements: python3

Program uses following moduli:

random, os, subprocess, copy, shutil, read_write_i

Input file - solution_structures.in with following parameters:

Description of the parameters see in example of the input file

Output data:

Structure files of solid solutions in rv_at format with corresponding visualization files in xyz format.
They can contain vacancies if structure type 'vacancy' was choosen.

Program can be used as importable or executable

\\newpage

'''

#====================================================================================
#================================= FUNCTIONS BLOCK ==================================
#====================================================================================
# Function of a directory creation
def create_dir(directory):
    '''
    Function is intended for make of a directory
    Input parameters:
    directory  -  name of directory
    Output parameters:
    Directory with nama as the argument
    '''
    import os
    import shutil as S
    try: os.mkdir(directory)
    except OSError: ...
    S.rmtree(directory)
    os.mkdir(directory)

# Function of building of the input file for md0_vol_defect and md0_vol_defect_small programs
def md0_vol_defect_in(pot_names, n_sort, name_at_ideal, name_at_defect, n_cycle, c_bulk, time_step, level):
    '''
    Function is intended for building of input file for md0_vol_defect_small program for MD modeling
    Input parameters:
    pot_names                          - list of potentials' names using in molecular dynamics (MD)
    n_sort                             - number of distinct elements
    name_at_ideal                      - name of the file with ideal defectless lattice
    name_at_defect                     - name of the file with defect lattice
    n_cycle_vol_defect                 - number of MD cycles
    new_cbulk                          - bulk modulus corresponding to the solute concentrations
    t_step                             - time step of MD (femtoseconds)
    level                              - how many times move to upper directory from calculation directory is needed to approach
                                         directory with potentials
    Output data:
    md0_vol_defect.in file for md0_vol_defect_small program
    '''
    l = ["---------------------------------\n",
         "5000                  : n_cylce_relax\n",
         "0.                    : pressure_ref_gpa\n",
         "170.                  : c_bulk ! in GPa (for refit pressure by Berendsen method)\n",
         "'dataw/id_r.at'    :(-) name_at_ideal\n",
         "'dataw/v1.at'   :(-) name_at_defect\n",
         "'dataw/id_r.at'  :(+) name_at_ideal_r\n",
         "'dataw/v1_r.at' :(+) name_at_defect_r\n",
         "3.                    : h (fs)\n",
         "1.                    : force_max (in meV/A)\n",
         "1.D-06                : pressure_max (in GPa)\n",
         "1.                  : dr_pot   ! Angstrom\n",
         "-5.31D0                 : e_at2_ideal ! if < 0 then calc h_sol_2in1\n", 
         "\n",
         "open(1,file='md0_vol_defect.in')\n",
         "read(1,*) name_pot\n",
         "read(1,*) n_cycle_relax\n", 
         "read(1,*) pressure_ref_gpa ! pressure in GPa\n",
         "read(1,*) c_bulk ! in GPa (for refit pressure by Berendsen method)\n",
         "read(1,*) name_at_ideal\n",
         "read(1,*) name_at_defect\n",
         "read(1,*) name_at_ideal_r\n",
         "read(1,*) name_at_defect_r\n",
         "read(1,*) h ! fs\n",
         "read(1,*) force_max\n",
         "read(1,*) pressure_max ! in GPa\n",
         "read(1,*) dr_pot   ! Angstrom\n",
         "read(1,*) e_at2_ideal ! for h_sol_2in1 if < 0\n",
         "close(1)\n"]
    pot_list = []
    for j in pot_names:
        pot_list.append("'"+level*"../"+"pot/"+j+"'          : \n")
    new_defect_list = [str(n_sort)+"                         : n_sort\n"]+pot_list+l
    new_defect_list[5] = n_cycle+"                  : n_cylce_relax\n"
    new_defect_list[7] = c_bulk+"                  : c_bulk ! in GPa (for refit pressure by Berendsen method)\n"
    new_defect_list[8] = "'dataw/"+name_at_ideal+"'          :(-) name_at_ideal\n"
    new_defect_list[9] = "'dataw/"+name_at_defect+"'         :(-) name_at_defect\n"
    new_defect_list[10] = "'dataw/"+name_at_ideal+"_r'       :(-) name_at_ideal_r\n"
    new_defect_list[11] = "'dataw/"+name_at_defect+"_r'      :(-) name_at_defect_r\n"
    new_defect_list[12] = time_step+"                    : h (fs)\n"
    f = open('md0_vol_defect.in', 'w')
    f.writelines(new_defect_list)
    f.close()

# Function of building of the scrpt for launching of md0_vol_defect and md0_vol_defect_small programs
def run_script(cluster, list_of_nodes):
    '''
    Function is intended for building of the script for launch of the md0_vol_defect_small program for MD modeling
    Input parameters:
    cluster                            - name of the computational cluster                           
    list_of_nodes                      - list of nodes' names on choosen computational cluster
    Output data:
    File of launch script of the md0_vol_defect_small program
    '''
    import random
    run_script_nsmn = [
        "#PBS -N md0_vol_def\n",
        "#PBS -l nodes=1:ppn=1\n",
        "#PBS -q q32gb\n",
        "#PBS -r n\n",
        "#PBS -o out\n",
        "#PBS -e err\n\n",

        "cd $PBS_O_WORKDIR\n\n",

        "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/build/gcc/4.6.2/lib:/home/build/gcc/4.6.2/lib64:/home/build/mpi/mpich2141/lib/trace_r\n\n",

        "/home/lipnic/bin/md0_vol_defect_small"
        ]

    run_script_knit = [
        "#PBS -N md0_vol_def\n",
        "#PBS -l nodes=node01:ppn=10+node02:ppn=9\n",
        "#PBS -q batch\n",
        "#PBS -r n\n",
        "#PBS -j oe\n",
        "#PBS -l walltime=99999999:00:00\n\n"
        "cd $PBS_O_WORKDIR\n\n",
        "/home/lipnickiy_a/bin/md0_vol_defect_small"
        ]
    if cluster == 'nsmn':
        cur_run_script = run_script_nsmn[:]
    elif cluster == 'knit':
        node = random.choice(list_of_nodes) 
        cur_run_script = run_script_knit[:]
        cur_run_script[1] = "#PBS -l nodes="+node+":ppn=1\n"
    f = open('run', 'w')
    f.writelines(cur_run_script)          
    f.close()

# Functions of matrix atom replacement by solute atoms
def insert_solute_atoms(num_el = 0, atomic_blocks = (), atomic_block_names = (), atomic_numbers = [], sol_at_concs = [], sol_at_sort_ats = [], sol_at_masses = []):
    '''
    This function is intended for random replacement of atoms in a pure element structure by solute atoms.
    Input parameters:
    num_el                             - number of distinct elements
    atomic_blocks                      - tuple of lists, length of list = number of atoms in block,
    [num_at_r[i], i_sort_at[i], mass_at[i], r_at[i], v_at[i]] (order number, sort, mass, 
    cartesian coordinates, velocity) - element of list
    atomic_block_names                 - tuple of block names
    atomic_numbers                     - list of atomic order numbers
    sol_at_concs                       - list of solute atoms concentrations
    sol_at_sort_ats                    - list of solute atoms sorts in final structure (e.g. for Ti-V-Fe system it will be [1,2,3])
    sol_at_masses                      - list of solute atoms masses
    Output data:
    One list, that is sum of all atomic blocks, with initial atoms replaced by solute elements
    '''
    import random
    new_structure = []
    blocks = {}
    all_solute_atoms = {}
    for i in range(len(atomic_blocks)):
        all_solute_atoms[i] = []
    for el in range(1, num_el):
        solute_concentration = sol_at_concs[el-1]
        new_i_sort_at = sol_at_sort_ats[el-1]
        solute_atom_mass = sol_at_masses[el-1]
        for l in range(len(atomic_blocks)):
            # Set solute concentrations in all blocks to choosen values
            cur_list = atomic_numbers[l] 
            n_at = len(atomic_numbers[l])
            n_subst = solute_concentration*n_at
            i_at = []
            for i in range(int(n_subst)):
                j = random.randint(0,n_at-1)              
                while (j in i_at) or (j in all_solute_atoms[l]): j = random.randint(0,n_at-1)              
                i_at.append(j)
            
            for i in range(n_at):
                if(i in i_at): 
                    num = atomic_numbers[l][i]  
                    for at in atomic_blocks[l]: 
                        if num == at[0]: 
                            at[1] = new_i_sort_at
                            at[2] = solute_atom_mass
            print('Atomic block '+str(l)+' len i_at = ', len(i_at))
            all_solute_atoms[l] += i_at

    # Build of one block as sum of initial blocks
    for l in range(len(atomic_blocks)):
        blocks[atomic_block_names[l]] = atomic_blocks[l]
        new_structure += atomic_blocks[l]
    insert_solute_atoms.blocks = blocks
    return new_structure

def selected_insert_solute_atoms(num_el = 0, at_names = [], atomic_blocks = (), atomic_block_names = (), atomic_numbers = (), sol_at_concs = [], sol_at_sort_ats = [], sol_at_masses = [], criterion = '', max_dist_cent_sol = 0, repl_at_lists= [], dict_main = {}):
    '''
    This function is intended for random replacement of atoms (choosen according to specific criterion) 
    in a pure element structure by solute atoms.
    Input parameters:
    num_el                             - number of distinct elements
    at_names                           - list of atom names corresponding to lists that are shown below
    atomic_blocks                      - tuple of lists, length of list = number of atoms in block,
    [num_at_r[i], i_sort_at[i], mass_at[i], r_at[i], v_at[i]] (order number, sort, mass, 
    cartesian coordinates, velocity) - element of list
    atomic_block_names                 - tuple of block names
    atomic_numbers                     - list of atomic order numbers
    sol_at_concs                       - list of solute atoms concentrations
    sol_at_sort_ats                    - list of solute atoms sorts in final structure (e.g. for Ti-V-Fe system it will be [1,2,3])
    sol_at_masses                      - list of solute atoms masses
    criterion                          - type of selection of the initial atoms for replacement 
                                         (current variants - 'distance_from_center' (center of initial lattice) - 
                                                              only atoms with distance from center less than choosen value
                                                              can be replaced by solute
                                                             'all' - all atoms can be replaced
                                                             'manual_replacement' - only choosen atoms will be replaced)
    max_dist_cent_sol                  - maximal distance between center of initial structure and initial atom 
                                         (used if criterion - 'distance_from_center')
    repl_at_lists                      - list of lists for manual replacement of atoms (if 'manual_replacement' is choosen) 
                                         (each list corresponds to element from at_names[1:] list, 
                                         at_names[0] - initial (matrix) element)
    dict_main                          - rv_at dictionary of initial structure
    Output data:
    One list, that is sum of all atomic blocks, with initial atoms replaced by solute elements
    (with replacement made by certain criterion)
    '''

    # Import of required moduli
    import copy
    import random

    # Dictionary of initial structure in rv_at format
    d = copy.deepcopy(dict_main)

    # Main part of the function
    new_structure = []
    blocks = {}
    at_list_dict = {}
    all_solute_atoms = {}
    for i in range(len(atomic_blocks)):
        all_solute_atoms[i] = []

    # Determination of the lattice center
    x_list = [a[0] for a in d['r_at']]
    y_list = [a[1] for a in d['r_at']]
    z_list = [a[2] for a in d['r_at']]
    x_center = (max(x_list)+min(x_list))/2
    y_center = (max(x_list)+min(x_list))/2
    z_center = (max(x_list)+min(x_list))/2

    # Subsequent replacement of initial atoms by solute elements in accordance to its concentrations in each block of atomic_blocks
    for l in range(len(atomic_blocks)):
        # Set solute concentrations in all blocks to choosen values
        cur_list = atomic_numbers[l] 
        n_at = len(atomic_numbers[l])
        block_name = atomic_block_names[l]
        at_list = []

        # Set of the criterion of initial atom choise for their random replacemen by solute elements
        if criterion == 'distance_from_center':
            if max_dist_cent_sol == 0: raise RuntimeError('Input maximum distance between the center of cell and solute atom!!!')
            for i in range(n_at):                
                r_at_cur = d['r_at'][i]
                dist_cur = ( (r_at_cur[0]-x_center)**2 + (r_at_cur[1]-y_center)**2 + (r_at_cur[2]-z_center)**2 )**0.5
                if dist_cur > max_dist_cent_sol:
                    continue 
                else: at_list.append(i)
            at_list_final = at_list

        elif criterion == 'all':
            at_list = list(range(n_at))
            at_list_final = at_list

        # Cycle by solute elemetns
        for el in range(1, num_el):
            # Set of parameters for current element
            solute_concentration = sol_at_concs[el-1]
            new_i_sort_at = sol_at_sort_ats[el-1]
            solute_atom_mass = sol_at_masses[el-1]
            at_name = at_names[el]
            n_subst = solute_concentration*n_at
            i_at = []

            # Replace of atoms by their numbers
            # List for replacement is different for each element
            if criterion == 'manual_replacement': 
                repl_at_list = repl_at_lists[el-1]
                if repl_at_list == []: 
                    raise RuntimeError('Input numbers of atoms for manual replacement!!!')
                for r in repl_at_list:
                    if r>(n_at-1): 
                        raise RuntimeError('Number cannot be greater than n_at-1!!!')
                at_list = [t-1 for t in repl_at_list]
                n_subst = len(repl_at_list)
                at_list_final = list(range(n_at))

            # Check for sites number in list of choosen atoms at_list
            if n_subst>len(at_list): 
                n_subst = len(at_list)
                print('Warning! All atoms on choosen sites from at_list will be replaced!')

            # Replacement of atoms
            for i in range(int(n_subst)):
                j = random.randint(0,len(at_list)-1)
                j1 = at_list[j]              
                while (j1 in i_at) or (j1 in all_solute_atoms[l]): 
                    j = random.randint(0,len(at_list)-1)
                    j1 = at_list[j]              
                i_at.append(j1)

            for i in range(n_at):
                if(i in i_at): 
                    num = atomic_numbers[l][i]  
                    for at in atomic_blocks[l]: 
                        if num == at[0]: 
                            at[1] = new_i_sort_at
                            at[2] = solute_atom_mass
            print('Atomic block '+block_name+' len i_at = ', len(i_at))
            all_solute_atoms[l] += i_at
            print('True concentration of element '+at_name+' in block '+block_name+' = '+str(len(i_at)/n_at))

    # Preparing of the output data
    for l in range(len(atomic_blocks)):
        blocks[atomic_block_names[l]] = atomic_blocks[l]
        new_structure += atomic_blocks[l]
    selected_insert_solute_atoms.blocks = blocks
    selected_insert_solute_atoms.at_list_final = at_list_final
    return new_structure

# Function of the vacancy introduction into solid solution on matrix atom sites
def selected_insert_vacancies(num_el=0, composition_name = '', version = None, at_names = [], matr_at_sort_at = 0, sol_at_sort_ats = [], sol_at_masses = [], criterion = '', max_dist_vac_sol = 0, repl_at_vac_list = [],
                              pot_names = [], n_sort = 0, n_cycle_vol_defect = 0, new_cbulk = 0, t_step = 0, level=0, md0_vol_defect_need = 'False',
                              cluster = '', list_of_nodes = []):
    '''
    This function is intended for replacement of matrix atoms by vacancies in solid solution. For now
    it subsequently replace only one matrix atom (over all matrix atoms), while other remain in lattice. 
    Input parameters:
    num_el                             - number of distinct elements
    composition_name                   - part of name of the structure, containing symbols of elements and their concentrations
    version                            - number of variant of the structure with equal composition but different atomic arrangement
    at_names                           - list of atom names corresponding to lists that are shown below
    matr_at_sort_at                    - sort of matrix (initial) atoms
    sol_at_sort_ats                    - list of solute atoms sorts in final structure (e.g. for Ti-V-Fe system it will be [1,2,3])
    sol_at_masses                      - list of solute atoms masses
    criterion                          - type of selection of the matrix atoms for replacement 
                                         (current variants - will be added later)
    max_dist_vac_sol                   - maximal distance between vacancy and any solute atom 
    repl_at_vac_list                   - list for manual replacement of matrix atoms by vacancy
                                         (each list corresponds to element from at_names[1:] list, 
                                         at_names[0] - initial (matrix) element)
    md0_vol_defect_need                - if 'True' then MD files will be built (all below parameters are relevant if this parameter
                                         is 'True')
    pot_names                          - list of potentials' names using in molecular dynamics (MD)
    n_sort                             - number of distinct elements
    n_cycle_vol_defect                 - number of MD cycles
    new_cbulk                          - bulk modulus corresponding to the solute concentrations
    t_step                             - time step of MD (femtoseconds)
    level                              - how many times move to upper directory from calculation directory is needed to approach
                                         directory with potentials
    cluster                            - name of the computational cluster                           
    list_of_nodes                      - list of nodes' names on choosen computational cluster
    Output data:
    Folder with composition_name, in which folders with names like '1', '2', etc. 
    (number of variant of structure with certain composition). These folders contain folders with names
    like '<number of replaced atom>_<minimal distance from any solute atom>'. Each of these folders contains
    files need for MD and visualized structures (in xyz format) - ideal and with vacancy. For clarity, on vacancy
    site nitrogen atom is placed.
    '''
    import solution_structures
    from read_write_i import ReadWrite as RW
    import shutil as S

    # Read of ideal structure
    e3 = RW()
    e3.r_rv_at('id_'+composition_name+'_'+str(version))
    d3 = e3.r_rv_at_dict

    # Determination of all matrix sites in solid solution
    matr_at_list = []
    sol_at_list = []
    sol_at_dict = {}
    for i in range(len(sol_at_sort_ats)):
        sol_at_dict[at_names[i+1]] = []

    for i in range(len(d3['r_at'])):
        if d3['i_sort_at'][i] == matr_at_sort_at: 
            matr_at_list.append(i)
        else:
            sol_at_list.append(i)
            for j in range(len(sol_at_sort_ats)):
                if d3['i_sort_at'][i] == sol_at_sort_ats[j]:  
                    sol_at_dict[at_names[j+1]].append(i)
    for key in sol_at_dict:
        print ('Number of atoms '+key+' = '+str(len(sol_at_dict[key])))

    # Creation of vacancies on each matrix site with saving of the derived structure
    for i in matr_at_list:
        d4 = copy.deepcopy(d3)
        e4 = RW()
        d4['n_at'] = d3['n_at']-1
        d4['mark_at'] = [['F' for f in range(10)] for f1 in range(d4['n_at'])]
        d4['mark_green'] = ['T' for f2 in range(d4['n_at'])]
        d4['num_at_r'] = [1+f3 for f3 in range(d4['n_at'])]
        d4['i_sort_at'] = []
        d4['mass_at'] = []
        d4['r_at'] = []
        d4['v_at'] = []                
        for j in range(len(d3['r_at'])):
            if j == i: continue
            d4['i_sort_at'].append(d3['i_sort_at'][j])
            d4['mass_at'].append(d3['mass_at'][j])                    
            d4['r_at'].append(d3['r_at'][j])
            d4['v_at'].append(d3['v_at'][j])

        # Determination of distance between vacancy and the closest solute atom
        distance_dict = {}
        for k in sol_at_list:
            distance = ((d3['r_at'][i][0]-d3['r_at'][k][0])**2+
                        (d3['r_at'][i][1]-d3['r_at'][k][1])**2+
                        (d3['r_at'][i][2]-d3['r_at'][k][2])**2)**0.5
            distance_dict[distance] = k
        min_dist = '{0:.3f}'.format(min(list(distance_dict.keys())))
        min_dist_float = float(min_dist)
        if min_dist_float > max_dist_vac_sol: continue
        e4.w_rv_at('vac_'+str(i+1)+'_'+str(version), d4)

        # Creation of folder for each variant of vacancy
        solution_structures.create_dir(composition_name+'/'+str(n)+'/'+str(i+1)+'_'+str(min_dist))
        solution_structures.create_dir(composition_name+'/'+str(n)+'/'+str(i+1)+'_'+str(min_dist)+'/dataw')

        if md0_vol_defect_need == 'True':
            # Build of input file and launch script for program md0_vol_defect_small
            name_id = 'id_'+composition_name+'_'+str(version)
            name_def = 'vac_'+str(i+1)+'_'+str(version)
            solution_structures.md0_vol_defect_in(pot_names, n_sort, name_id, name_def, n_cycle_vol_defect, str(new_cbulk), t_step, level)
            solution_structures.run_script(cluster, list_of_nodes)

        # Visualization of the structure in xyz format
        f1 = open('vac_'+str(i+1)+'_'+str(version)+'.xyz', 'w')
        f1.write(str(d4['n_at']+1)+'\n')
        f1.write('# Crystal structure\n')
        f1.write('{0:5s}'.format('N')+
                     '{0:12.6f}'.format(d3['r_at'][i][0])+
                     '{0:12.6f}'.format(d3['r_at'][i][1])+
                     '{0:12.6f}'.format(d3['r_at'][i][2])+'\n')
        for at in range(len(d4['r_at'])):
            if d4['i_sort_at'][at] == matr_at_sort_at: 
                symb = at_names[0]
            else: 
                for el in range(1, num_el):
                    s_at = sol_at_sort_ats[el-1]
                    if d4['i_sort_at'][at] == s_at:
                        symb = at_names[el]
            f1.write('{0:5s}'.format(symb)+
                     '{0:12.6f}'.format(d4['r_at'][at][0])+
                     '{0:12.6f}'.format(d4['r_at'][at][1])+
                     '{0:12.6f}'.format(d4['r_at'][at][2])+'\n')
        f1.close()

        # Move of files into corresponding directories
        S.copy('id_'+composition_name+'_'+str(version), composition_name+'/'+str(n)+'/'+str(i+1)+'_'+str(min_dist)+'/dataw')
        S.move('vac_'+str(i+1)+'_'+str(version),        composition_name+'/'+str(n)+'/'+str(i+1)+'_'+str(min_dist)+'/dataw')
        S.move('vac_'+str(i+1)+'_'+str(version)+'.xyz', composition_name+'/'+str(n)+'/'+str(i+1)+'_'+str(min_dist))
        S.move('md0_vol_defect.in',                     composition_name+'/'+str(n)+'/'+str(i+1)+'_'+str(min_dist))
        S.move('run',                                   composition_name+'/'+str(n)+'/'+str(i+1)+'_'+str(min_dist))

# Function of lattice scaling according to the solute atom concentrations
def lattice_scaling(num_el = 0, sol_at_concs = [], sol_at_acells = [], matr_at_acell = 0, latt_dict = {}):
    '''
    Function is intended for scaling of the lattice parameters according to the solute atoms concentrations
    Input parameters:
    num_el                             - number of distinct elements
    sol_at_concs                       - list of solute atoms concentrations
    sol_at_acells                      - list of lattice parameters of structures with pure solute atoms at their ground states
    matr_at_acell                      - lattice parameter of structure with pure matrix atom at its ground state
    latt_dict                          - rv_at dictionary of the solid solution structure
    Output data:
    rv_at dictionary with scaled lattice parameters, sizes and atomic coordinates
    '''
    import copy
    d = copy.deepcopy(latt_dict)
    full_solute_concentration = 0
    full_solute_acell = 0
    for el in range(1, num_el):
        full_solute_concentration += sol_at_concs[el-1]
        full_solute_acell += sol_at_concs[el-1]*sol_at_acells[el-1]
    new_acell = matr_at_acell*(1-full_solute_concentration) + full_solute_acell
    coef = new_acell/d['a_lattice3'][0]
    d['size'] = [i*coef for i in d2['size']]
    d['a_lattice3'] = [i*coef for i in d2['a_lattice3']]
    d['r_at'] = [[i[0]*coef, i[1]*coef, i[2]*coef] for i in d2['r_at']]
    return d

# Function of bulk modulus scaling according to the solute atom concentrations
def cbulk_scaling(num_el = 0, sol_at_concs = [], sol_at_cbulks = [], matr_at_cbulk = 0):
    '''
    Function is intended for scaling of the bulk modulus according to the solute atoms concentrations
    Input parameters:
    num_el                             - number of distinct elements
    sol_at_concs                       - list of solute atoms concentrations
    sol_at_cbulks                      - list of bulk moduli of structures with pure solute atoms at their ground states
    matr_at_cbulk                      - bulk modulus of the structure with pure matrix atom at its ground state
    Output data:
    Bulk modulus of solid solution structure
    '''
    full_solute_concentration = 0
    full_solute_cbulk = 0
    for el in range(1, num_el):
        full_solute_concentration += sol_at_concs[el-1]
        full_solute_cbulk += sol_at_concs[el-1]*sol_at_cbulks[el-1]
    new_cbulk = matr_at_cbulk*(1-full_solute_concentration) + full_solute_cbulk
    return new_cbulk

# Function of structure visualization in xyz format
def xyz_format(num_el = 0, atomic_block = [], name_block = '', name_composition = '', at_names = [], matr_at_sort_at = 0, sol_at_sort_ats = []):
    '''
    Function is intended for visualization of the structure in xyz format
    Input parameters:
    num_el                             - number of distinct elements
    atomic_block                       - list with length = number of atoms in block,
    [num_at_r[i], i_sort_at[i], mass_at[i], r_at[i], v_at[i]] (order number, sort, mass, 
    cartesian coordinates, velocity) - element of list
    name_block                         - name of the atomic_block
    name_composition                   - part of name of the structure, containing symbols of elements and their concentrations
    at_names                           - list of atom names 
    matr_at_sort_at                    - sort of matrix (initial) atoms
    sol_at_sort_ats                    - list of solute atoms sorts in final structure (e.g. for Ti-V-Fe system it will be [1,2,3])
    Output data:
    File in xyz format corresponding to structure in atomic_block
    '''
    if matr_at_sort_at == 0: raise RuntimeError('Atomic sort for matrix atoms was not set!!!')
    symbs = []
    exec(name_block+'_count_'+at_names[0]+'= 0')
    symbs.append(at_names[0])
    for j in range(1, num_el):
        s = at_names[j]
        exec(name_block+'_count_'+s+' = 0')
        symbs.append(s)
    f1 = open(name_block+'_'+name_composition+'.xyz', 'w')
    f1.write(str(len(atomic_block))+'\n')
    f1.write('# Crystal structure\n')
    for at in atomic_block:
        if at[1] == matr_at_sort_at: 
            symb = at_names[0]
            exec(name_block+'_count_'+symb+' += 1')
        else:
            for el in range(1, num_el):
                s_at = sol_at_sort_ats[el-1]
                if at[1] == s_at:
                    symb = at_names[el]
                    exec(name_block+'_count_'+symb+' += 1')
        f1.write('{0:5s}'.format(symb)+'{0:12.6f}'.format(at[3][0])+'{0:12.6f}'.format(at[3][1])+'{0:12.6f}'.format(at[3][2])+'\n')

    for element in symbs:
        print(name_block+'_count_'+str(element)+' = ', eval(name_block+'_count_'+str(element)))
    f1.close()

# Function of building of the string, containing symbols of elements and its concentrations in solid solution
def name_composition(num_el = 0, at_names = [], sol_at_concs = [], version = None):
    '''
    Function is intended for building of part of the structures' names, containing symbols of the elements and their concentrations
    Input parameters:
    num_el                             - number of distinct elements
    at_names                           - list of atom names corresponding to lists that are shown below
    sol_at_concs                       - list of solute atoms concentrations
    version                            - number of variant of the structure with equal composition but different atomic arrangement
    Output data:
    E.g. for V-0.2atTi-0.4at.Fe it will be 'V_Ti_0.2_Fe_0.4' (if version = None) or 'V_Ti_0.2_Fe_0.4_3' (if version = 3)
    '''
    name = at_names[0]+'_'
    for el in range(1, num_el):
        name += at_names[el]+'_'+str(sol_at_concs[el-1])+'_'
    if version != None: name += str(version)
    else: name=name[:-1]
    return name

# Function of writing of solute atoms into origin pure structure
def write_replaced_atoms(changes = [], dict_aux = {}, dict_main = {}):
    '''
    Function is intended for writing of the changed structure given in format of atomic_block
    (which is list of following elements[num_at_r[i], i_sort_at[i], mass_at[i], r_at[i], v_at[i]] (order number, sort, mass, 
    cartesian coordinates, velocity)) in rv_at dictionary
    Input parameters:
    changes    -   sum of all atomic blocks
    dict_aux   -   rv_at dictionary of the initial structure
    dict_main  -   rv_at dictionary for writing changes
    Output data:
    rv_at dictionary of solid solution
    '''
    count = 1
    dict_main['num_at_r'] = []
    dict_main['i_sort_at'] = []
    dict_main['mass_at'] = []
    dict_main['r_at'] = []
    dict_main['v_at'] = []
    for i in range(len(dict_aux['r_at'])):
        dict_main['num_at_r'].append(count)
        dict_main['i_sort_at'].append(changes[i][1])
        dict_main['mass_at'].append(changes[i][2])
        dict_main['r_at'].append(changes[i][3])
        dict_main['v_at'].append(changes[i][4])
        count += 1
    print('Total number of atoms = ', count-1)
    return dict_main

# Function of atomic block formation from full structure (written in rv_at format) dictionary
def extract_block(dict_main = {}, criterion = 'True'):
    '''
    Function is intended for building of the atomic_block from rv_at dictionary of initial structure
    by certain criterion
    Input parameters:
    dict_main  -  rv_at dictionary of the initial structure
    criterion  -  condition for separating of the initial structure into atomic blocks (e.g. z (coordination) > 0)
    Output data:
    Atomic block satisfying to the criterion.
    Atomic_block is list of following elements[num_at_r[i], i_sort_at[i], mass_at[i], r_at[i], v_at[i]] (order number, sort, mass, 
    cartesian coordinates, velocity)
    '''
    block = []
    block_numbers = []
    for i in range(len(d1['r_at'])):
        if eval(criterion): 
            block.append([d1['num_at_r'][i], 
                             d1['i_sort_at'][i], 
                             d1['mass_at'][i], 
                             d1['r_at'][i], 
                             d1['v_at'][i]])
            block_numbers.append(d1['num_at_r'][i])
    return block, block_numbers

if __name__ == '__main__':

    import random
    import subprocess as SP
    import copy
    import os
    import shutil as S
    from read_write_i import ReadWrite as RW

    #====================================================================================
    #================================ READING INPUT DATA ================================
    #====================================================================================
    
    # Determination of the current directory
    home_dir = os.getcwd()

    # Initial values of some parameters
    replace_atom_list = []
    atom_names = []
    elements = {}
    pot_names = []
    n_sort = 0

    # Read input file
    f = open('solution_structures.in')
    l = f.readlines()
    f.close()
    f = open('solution_structures.in')
    for line in f:
        line1 = line.split()
        if len(line) == 1 or line[0] == '#': continue

        # Common parameters for all variants of type_of_structure
        elif 'type_of_structure' == line1[0]: type_of_structure = line1[1]
        elif 'number_of_elements' == line1[0]:
            number_of_elements = int(line1[1])
        elif 'matrix_atom' == line1[0]: 
            matrix_atom_name = line1[1]
            matrix_atom_mass = line1[2]
            matrix_atom_acell = float(line1[3])
            matrix_atom_cbulk = float(line1[4])
            matrix_atom_sort_at = int(line1[5])
            atom_names.append(matrix_atom_name)
            elements[matrix_atom_sort_at] = matrix_atom_name
        elif 'switch_i_sort_at' == line1[0]: switch_i_sort_at = line1[1]
        elif 'max_distance_center_solute' == line1[0]: max_distance_center_solute = float(line1[1])
        elif 'number_of_solution_configs' == line1[0]: number_of_solution_configs = int(line1[1])
        elif 'path_to_initial_file' == line1[0]: path_to_initial_file = line1[1]

        # Parameters for type_of_structure = 'melting'
        elif 'interface_position_z' == line1[0]: interface_position_z = float(line1[1]) 

        # Computational cluster parameters
        elif 'cluster' == line1[0]: cluster = line1[1]
        elif 'list_of_nodes' == line1[0]: list_of_nodes = line1[1].split(',')

        # Molecular dynamics calculations parameters
        elif 'md0_vol_defect_need' == line1[0]: md0_vol_defect_need = line1[1]
        elif 'pot1' == line1[0]: 
            pot1 = line1[1]
            pot_names.append(pot1)
            n_sort += 1
        elif 'pot2' == line1[0]: 
            pot2 = line1[1]
            pot_names.append(pot2)
            n_sort += 1
        elif 'pot12' == line1[0]: 
            pot12 = line1[1]
            pot_names.append(pot12)
        elif 'pot12' == line1[0]: pot12 = line1[1]
        elif 'max_distance_vacancy_solute' == line1[0]: max_distance_vacancy_solute = float(line1[1])
        elif 'n_cycle_vol_defect' == line1[0]: n_cycle_vol_defect = line1[1]
        elif 't_step' == line1[0]: t_step = line1[1]
        elif 'level' == line1[0]: level = int(line1[1])

    f.close()

    # Check for solute atoms presence for formation of the solution
    if number_of_elements == 1: raise RuntimeError('Your structure contains only one element!!!')

    # Read of input data for solute elements
    solute_atom_concentrations = []
    solute_atom_masses = []
    solute_atom_sort_ats = []
    solute_atom_acells = []
    solute_atom_cbulks = []
    replace_atom_lists = []
    for i in range(1, number_of_elements):
        for line in l:
            line1 = line.split()
            if 'solute_atom_'+str(i) in line: 
                print('Sabaton')
                exec('solute_atom_name_'+str(i)+' = line1[1]')
                exec('solute_atom_mass_'+str(i)+' = line1[2]')
                exec('solute_atom_acell_'+str(i)+' = float(line1[3])')
                exec('solute_atom_cbulk_'+str(i)+' = float(line1[4])')
                exec('solute_atom_sort_at_'+str(i)+' = int(line1[5])')
                exec('solute_atom_concentration_'+str(i)+' = float(line1[6])')
                exec('solute_atom_concentrations.append(solute_atom_concentration_'+str(i)+')' )
                exec('solute_atom_masses.append(solute_atom_mass_'+str(i)+')' )
                exec('solute_atom_sort_ats.append(solute_atom_sort_at_'+str(i)+')' )
                exec('solute_atom_acells.append(solute_atom_acell_'+str(i)+')' )
                exec('solute_atom_cbulks.append(solute_atom_cbulk_'+str(i)+')' )
                exec('atom_names.append(solute_atom_name_'+str(i)+')' )
                exec('elements[solute_atom_sort_at_'+str(i)+'] = solute_atom_name_'+str(i))
            elif 'replace_atom_list_'+str(i) in line:
                exec('replace_atom_list_'+str(j)+' = [int(j)-1 for j in line1[1:]]')

    # Read of initial structure file into dictionary
    e1 = RW()
    e1.r_rv_at(path_to_initial_file)
    d1 = e1.r_rv_at_dict
    if switch_i_sort_at == 'True': 
        d1['i_sort_at']=[matrix_atom_sort_at for i in range(len(d1['i_sort_at']))]


    #====================================================================================
    #= BUILD OF THE STRUCTURE FOR DETERMINATION OF MELTING POINT OF THE SOLID SOLUTION ==
    #====================================================================================
    if type_of_structure == 'melting':
        # Split of full number of atoms on two parts relative to the interface position.
        # Each part is list of lists, containing following information about atom:
        # [num_at_r[i], i_sort_at[i], mass_at[i], r_at[i], v_at[i]] (order number, sort, mass, 
        # cartesian coordinates, velocity)
        # Parameters n_at, n_mark_at, size, a_lattice3, mark_at, mark_green are common for both.
        main_top_part, main_top_part_numbers = extract_block(dict_main = d1, 
                                                             criterion = 'd1["r_at"][i][2]>0')
        main_low_part, main_low_part_numbers = extract_block(dict_main = d1, 
                                                             criterion = 'd1["r_at"][i][2]<=0')

        # Random choise of atoms for replacement by solute atoms with subsequent writing of 
        # the structure into file in rv_at format
        for n in range(1,number_of_solution_configs+1):
            # Preparing of current lists and dictionaries
            d2 = copy.deepcopy(d1)
            top_part = copy.deepcopy(main_top_part)
            top_part_numbers = copy.deepcopy(main_top_part_numbers)
            low_part = copy.deepcopy(main_low_part)
            low_part_numbers = copy.deepcopy(main_low_part_numbers)

            # Replacement of matrix atoms
            changed_structure = selected_insert_solute_atoms(num_el = number_of_elements,
                                                             at_names = atom_names,
                                                             atomic_blocks = (top_part, low_part), 
                                                             atomic_block_names = ('top_part', 'low_part'), 
                                                             atomic_numbers = (top_part_numbers, low_part_numbers), 
                                                             sol_at_concs = solute_atom_concentrations, 
                                                             sol_at_sort_ats = solute_atom_sort_ats, 
                                                             sol_at_masses = solute_atom_masses, 
                                                             criterion = 'all', 
                                                             max_dist_cent_sol = max_distance_center_solute, 
                                                             repl_at_lists = [], 
                                                             dict_main = d2)

            # Writing of changes into rv_at dictionary
            d2 = write_replaced_atoms(changes = changed_structure, 
                                      dict_aux = d1, 
                                      dict_main = d2)

            # Lattice parameters scaling according to the solute concentrations
            d2 = lattice_scaling(num_el = number_of_elements, 
                                 sol_at_concs = solute_atom_concentrations, 
                                 sol_at_acells = solute_atom_acells, 
                                 matr_at_acell = matrix_atom_acell, 
                                 latt_dict = d2)

            # Writing rv_at file with solid solution
            name = name_composition(num_el = number_of_elements, 
                        at_names = atom_names, 
                        sol_at_concs = solute_atom_concentrations, 
                        version = n)
            e1.w_rv_at('at_melt_'+name, d2)

            # Visualization of all parts of the structure in xyz format
            blocks = selected_insert_solute_atoms.blocks
            for key in blocks:
                xyz_format(num_el = number_of_elements, 
                           atomic_block = blocks[key], 
                           name_block = key, 
                           name_composition = name, 
                           at_names = atom_names, 
                           matr_at_sort_at = matrix_atom_sort_at, 
                           sol_at_sort_ats = solute_atom_sort_ats)

    #====================================================================================
    #========================= VACANCY FORMATION IN SOLID SOLUTION ======================
    #====================================================================================
    elif type_of_structure == 'vacancy':

        # Build of the directory for choosen concentrations of solute elements
        name = name_composition(num_el = number_of_elements, 
                                at_names = atom_names, 
                                sol_at_concs = solute_atom_concentrations, 
                                version = None)
        create_dir(name)

        #======================== BUILD OF THE IDEAL SOLID SOLUTION =====================
        # Conversion of set of the atoms into list of lists, containing following information about atom:
        # [num_at_r[i], i_sort_at[i], mass_at[i], r_at[i], v_at[i]] (order number, sort, mass, 
        # cartesian coordinates, velocity)
        # Parameters n_at, n_mark_at, size, a_lattice3, mark_at, mark_green are common for both.
        main_top_part, main_top_part_numbers = extract_block(dict_main = d1, 
                                                             criterion = 'True')

        # Determination of the bulk modulus according to solute concentrations
        new_cbulk = cbulk_scaling(num_el = number_of_elements, 
                                  sol_at_concs = solute_atom_concentrations, 
                                  sol_at_cbulks = solute_atom_cbulks, 
                                  matr_at_cbulk = matrix_atom_cbulk)

        # Random choise of atoms for replacement by solute atoms with subsequent writing of 
        # the structure into file in rv_at format      
        for n in range(1,number_of_solution_configs+1):
            create_dir(name+'/'+str(n))
            # Preparing of current lists and dictionaries
            d2 = copy.deepcopy(d1)
            id_lat = copy.deepcopy(main_top_part)
            id_lat_numbers = copy.deepcopy(main_top_part_numbers)


            # Replacement of matrix atoms
            changed_structure = selected_insert_solute_atoms(num_el = number_of_elements,
                                                             at_names = atom_names,
                                                             atomic_blocks = (id_lat,), 
                                                             atomic_block_names = ('id',), 
                                                             atomic_numbers = (id_lat_numbers,), 
                                                             sol_at_concs = solute_atom_concentrations, 
                                                             sol_at_sort_ats = solute_atom_sort_ats, 
                                                             sol_at_masses = solute_atom_masses, 
                                                             criterion = 'distance_from_center', 
                                                             max_dist_cent_sol = max_distance_center_solute, 
                                                             repl_at_lists = [], 
                                                             dict_main = d2)
            # Writing of changes into rv_at dictionary
            d2 = write_replaced_atoms(changes = changed_structure, 
                                      dict_aux = d1, 
                                      dict_main = d2)

            # Lattice parameters scaling according to the solute concentrations
            d2 = lattice_scaling(num_el = number_of_elements, 
                                 sol_at_concs = solute_atom_concentrations, 
                                 sol_at_acells = solute_atom_acells, 
                                 matr_at_acell = matrix_atom_acell, 
                                 latt_dict = d2)
            e1.w_rv_at('id_'+name+'_'+str(n), d2)

            # Visualization of ideal solid solution in xyz format
            blocks = selected_insert_solute_atoms.blocks
            for key in blocks:
                xyz_format(num_el = number_of_elements, 
                           atomic_block = blocks[key], 
                           name_block = key, 
                           name_composition = name+'_'+str(n), 
                           at_names = atom_names, 
                           matr_at_sort_at = matrix_atom_sort_at, 
                           sol_at_sort_ats = solute_atom_sort_ats)

            # Copy of file with the ideal structure and its xyz representation into corresponding folder
            S.copy('id_'+name+'_'+str(n), name+'/'+str(n))
            S.copy('id_'+name+'_'+str(n)+'.xyz', name+'/'+str(n))
        
            # ============= INSERT OF VACANCIES INTO SOLID SOLUTION ====================
            selected_insert_vacancies(num_el = number_of_elements,
                                      composition_name = name, 
                                      version = n, 
                                      at_names = atom_names, 
                                      matr_at_sort_at = matrix_atom_sort_at, 
                                      sol_at_sort_ats = solute_atom_sort_ats, 
                                      sol_at_masses = solute_atom_masses, 
                                      criterion = '', 
                                      max_dist_vac_sol = max_distance_vacancy_solute, 
                                      repl_at_vac_list = [],
                                      pot_names = pot_names,
                                      n_sort = n_sort, 
                                      n_cycle_vol_defect = n_cycle_vol_defect, 
                                      new_cbulk = new_cbulk, 
                                      t_step = t_step, 
                                      level=level,
                                      md0_vol_defect_need = md0_vol_defect_need,
                                      cluster = cluster, 
                                      list_of_nodes = list_of_nodes)


