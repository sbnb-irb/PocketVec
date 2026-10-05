from subprocess import Popen, PIPE
import scipy.stats as ss
from tqdm import tqdm
import pandas as pd
import numpy as np
from Bio.PDB import *
import pickle
import shutil
import copy
import pybel
import sys
import os
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile


def prepare_pdb_moe(infile, outfile, logfile, moebatch = "/aloy/home/acomajuncosa/programs/MOE/moe2020/bin/moebatch"):
    """
    Prepare a PDB file for rDock docking. Requires MOE license & bin files.
    
    Args:
    
        infile (str): Path to input structure file (PDB format)
        outfile (str): Path to output structure file (MOL2 format)
        logfile (str): Path to LOG file. Check the logs for further details on protein preparation.
        moebatch (str): Path to your moebatch file. 
        
        
        

    """
    
    # Paths to moe functions
    moefunct = "/aloy/home/acomajuncosa/PocketVec_v2/GitLab_repo/code/utils/moefunctions.svl"
    
    # Prepare structures 
    command = moebatch + " -load " + moefunct + " -exec \"proteinprep['" + infile + "','" + outfile + "','" + logfile + "']\""
    
    process=Popen([command], stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = process.communicate()
    
    if "error" not in str(stderr).lower():
        sys.stderr.write("Preparation completed!")
        
    else:
        sys.stderr.write(str(stdout) + "\n\n")
        sys.stderr.write(str(stderr) + "\n\n")
        

def prepare_pdb_pdbfixer(infile, outfile, logfile):
    """
    Prepare a PDB file for rDock docking using the OpenSource software PDBFixer.
    
    Args:
    
        infile (str): Path to input structure file (PDB format)
        outfile (str): Path to output structure file (MOL2 format)
        logfile (str): Path to LOG file. 
    """

    try:
        fixer = PDBFixer(filename = infile)
        fixer.removeHeterogens(False)
        
        fixer.missingResidues = {}
        
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        
        fixer.findMissingAtoms()
        fixer.missingTerminals = {}
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        pdb = infile.split("/")[-1].split("_")[0]
        with open(logfile, "w") as log:
            log.write("Structure: " + pdb + "\n")
            log.write("·······················································" + "\n")
            
            log.write("Non-standard residues: " + str(len(fixer.nonstandardResidues)) + "\n")
            for res in fixer.nonstandardResidues:
                log.write(str(res[0]) + " -> " + str(res[1]) + "\n")
                
            log.write("·······················································" + "\n")
        
            log.write("Residues missing atoms: " + str(len(fixer.missingAtoms)) + "\n")
            for res in fixer.missingAtoms.keys():
                log.write(str(res) + " missing " + str(len(fixer.missingAtoms[res])) + " atoms \n")
                for atom in fixer.missingAtoms[res]:
                    log.write("\t" + str(atom) + "\n")
        
            log.write("·······················································" + "\n")

        pdb_file = outfile.replace(".mol2", ".pdb")
        PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'))
        mol = next(pybel.readfile("pdb", pdb_file))
        mol.write("mol2", outfile, overwrite=True)

        sys.stderr.write("Preparation completed")

    except Exception as e:
        sys.stderr.write(f"Error while preparing: {str(e)}")
        
        
def create_parameter_file(outfile, path_to_st, path_to_ctr, radius=str(12.0)):
    """
    Create structure parameter file.
    
    Args:
    
        outfile (str): Path to output structure parameter file (rDock-PRM format)
        path_to_st (str): Path to receptor (MOL2 format).
        path_to_ctr (str): Path to pocket centroid (SD format). 
        radius (float/str): Pocket radius (default 12.0).
        
    
    """
    
    text = """RBT_PARAMETER_FILE_V1.00
TITLE Change this if you want

RECEPTOR_FILE """ + path_to_st + """
RECEPTOR_FLEX 3.0

#################################################################
## CAVITY DEFINITION: REFERENCE LIGAND METHOD
#################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL """ + path_to_ctr + """
    RADIUS """ + radius + """
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
END_SECTION


SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION

#################################################
# WARNING!!! WE ARE DOING RIGID LIGAND DOCKING!!
#################################################
SECTION LIGAND
    TRANS_MODE FREE
    ROT_MODE FREE
    DIHEDRAL_MODE """ + "FIXED" + """
END_SECTION"""
    
    with open(outfile, "w") as f:
        f.write(text)
        
        
        
        
def environmental_variables(path_to_rDock):
    """
    Set environmental variables to run rDock.
    
    Args:
       
        path_to_rDock (str): Path to rDock installation directory. You may use the precompiled files we
        provide in code/utils/rDock_compiled or your local installed version of rDock.
    
    
    """
    
    if path_to_rDock[-1] == '/':
        path_to_rDock = path_to_rDock[:-1]
    
    os.environ["RBT_ROOT"] = path_to_rDock
    os.environ["LD_LIBRARY_PATH"] = path_to_rDock + "/lib"
    os.environ["PATH"] = os.environ["PATH"] + ":" + path_to_rDock + "/bin"
    
    
    

    
def create_cavity(outpath, path_to_log, path_to_cavity, path_to_rbcavity='rbcavity'):
    """
    Create cavity before docking with rDock. 3 files are generated:
        - log file
        - .grd file (docking cavity)
        - .as file
    
    Args:
        
        path_to_parameters (str): Path to structure parameters file
        path_to_log (str): Path to log
        path_to_cavity (str): Path to cavity (GRD extension)
        path_to_rbcavity (str): Path to rbcavity (if environment variables are set properly, 'rbcavity' should be enough)
    
    
    """
    
    # Path to parameters file
    path_to_parameters = os.path.join(outpath, 'st_parameters.prm')
    
    # Run rbcavity
    command = path_to_rbcavity + ' -was -d -r ' + path_to_parameters + " | tee " + path_to_log
    out = os.popen(command).read()
    sys.stderr.write(out)
    
    # Rename grid file
    os.rename(os.path.join(outpath, "st_parameters_cav1.grd"), os.path.join(outpath, "cavity.grd"))
    
    
    
    
    
    
def run_rDock(outpath, path_to_lib='../data/libs/TOP_128_rDock_LLM.sdf', nruns=25, seed=42):
    
    """
    Run rDock. Results are placed in the 'results.sd' file located in 'outpath'.
    
    Args:
        
        outpath (str): Path where all parameter files (2 x .prm) are located. All generated files will also be located in this path.
        path_to_lib (str): Path to the library of molecules to dock. By default, 128 LLMs are used.
        nruns (int/str): Number of docking runs for each molecule. By default, 25 runs are performed.
        seed (int/str): Seed to make docking calculations reproducible.
    
    
    """

    # Define variables
    results = os.path.join(outpath, 'results')  # SD file with results
    path_to_st_parameters = os.path.join(outpath, 'st_parameters.prm') # Parameters (st)
    path_to_rdock_parameters = os.path.join(outpath, 'dock.prm') # Parameters (dock)
    rdock_logs = os.path.join(outpath, 'rDock_log.log')
    nruns = str(nruns)
    seed = str(seed)

    command = 'rbdock -i ' + path_to_lib + ' -o ' + results + ' -r ' + path_to_st_parameters + ' -p ' + path_to_rdock_parameters + ' -n ' + str(nruns) + ' -s ' + seed + ' | tee ' + rdock_logs
    os.system(command)
    
    
    
    
def create_file_scores(path_to_scores, path_to_results):
    
    """
    Given a file with docking results (results.sd), create a file summarizing the minimum score for each docked molecule.
    
    Args:
        
        path_to_scores (str): Summary file with minimum docking scores.
        path_to_results (str): File with docking results (sd format)
    
    """
    
    
    with open(path_to_scores, "w") as out_score:
        
        # Read file
        file = pybel.readfile("sd", path_to_results)
        
        # Min affinity is set to an arbitrary large number
        aff = 1000000
        affinities = {}
        
        # Iterate over molecules and store the minimum score for each one
        for molecule in file:
            if molecule.title not in affinities:
                affinities[molecule.title] = aff
            value = float(molecule.data['SCORE.INTER'])
            if value < affinities[molecule.title]:
                affinities[molecule.title] = value

        # Print results
        for molecule in sorted(affinities):
            out_score.write(molecule + "\t" + str(affinities[molecule]) + "\n")
    

def rank_fp(fp_raw):
    
    """
    Given a fp with raw docking scores, convert it to ranks and deal with outlier molecules appropriately
    
    Args:
    
        fp_raw (list/np array): fp with raw docking scores
    
    """
    
    # Rank
    fp = ss.rankdata(fp_raw, method='max')
    # Increase ranking for positive scores
    for c in range(len(fp)):
        if fp_raw[c] > 0 and fp_raw[c] < 50:
            fp[c] = len(fp) + 1
        elif fp_raw[c] > 50 and fp_raw[c] < 100:
            fp[c] = len(fp) + 2
        elif fp_raw[c] > 100:
            fp[c] = len(fp) + 3
    return fp
    

def read_rDock_scores(path_to_scores):
    
    """
    Function to read tab-separated scores (molecule \t score)
    
    Args:
    
        path_to_scores (str): Path to file.
    
    """
    
    scores = {}
    with open(path_to_scores, "r") as f:
        for l in f:
            scores[l.split("\t")[0]] = float(l.split("\t")[1])
    return scores




def raw_fp(dict_scores, file_order):
    
    """
    Given a dict with the minimum docking score for each molecule and file with the order of molecules, create a fp of raw docking scores
    
    Args:
    
        dict_scores (dict): molecule-docking score
        file_order (str): path to a pickle file having a list with the ordered molecules
    
    """
    
    molecules = pickle.load(open(file_order, "rb"))
    return np.array([dict_scores[i] for i in molecules])


def select_first_model(path_in, path_out, is_mmcif):
    """
    
    Take only the first model of a PDB file (e.g NMR structure)
    
    Args:
       
        path_in (str): Path to input file
        path_out (str): Path to output file
        is_mmcif (boolean): Indicates if it is an mmCIF or not. 
    
    """
    if not is_mmcif:
        structure = PDBParser(QUIET=True).get_structure("st", path_in)
        io = PDBIO()

    else:
        structure = MMCIFParser(QUIET=True).get_structure("st", path_in)
        io = MMCIFIO()

    if len(structure) > 1:
        structure = structure[0]

    io.set_structure(structure)
    io.save(path_out)
    
def select_chain(path_in, path_out, chain, is_mmcif):
    """
    
    Take only the specified chain
    
    Args:
       
        path_in (str): Path to input file
        path_out (str): Path to output file
        chain (str): Chain ID
        is_mmcif (boolean): Indicates if it is an mmCIF or not. 
    
    """

    class ChainSelect(Select):
        def __init__(self, chain_id):
            self.chain_id = chain_id

        def accept_residue(self, residue):
            if residue.get_parent().id == self.chain_id:
                return 1
            else:
                return 0

    if not is_mmcif:
        structure = PDBParser(QUIET=True).get_structure("st", path_in)
        io = PDBIO()
        io.set_structure(structure)
        io.save(path_out, select_chain())

    else:
        structure = MMCIFParser(QUIET=True).get_structure("st", path_in)
        for residue in structure.get_residues():
            if len(residue.resname) > 3:
                residue.resname = residue.resname[:3]

        io = MMCIFIO()
        io.set_structure(structure)
        path_intermediate = os.path.join(path_out[:-4] + "_intermediate.cif")
        io.save(path_intermediate, ChainSelect(chain))


        structure = MMCIFParser(QUIET=True).get_structure("st", path_intermediate)
        for c in structure.get_chains():
            c.id = chain[0]
            break

        io = PDBIO()
        io.set_structure(structure)
        io.save(path_out, ChainSelect(chain[0]))
        
        os.remove(path_intermediate)
    
    
def remove_ligands(path_in, path_out, ligands):
    """
    
    Remove ligands from structure.
    
    Args:
       
        path_in (str): Path to input file
        path_out (str): Path to output file
        ligands (list/set): Ligands/residues to remove
        
    
    """

    structure = PDBParser(QUIET=True).get_structure("st", path_in)

    class remove_ligand(Select):
        def accept_residue(self, residue):
            if residue.get_resname() not in ligands:
                return 1
            else:
                return 0

    io = PDBIO()
    io.set_structure(structure)
    io.save(path_out, remove_ligand())
    
    
def remove_hydrogens(path_to_check_structure, path_in, path_out):
    """
    
    Remove hydrogens from residues
    
    Args:
       
        path_to_check_structure (str): Path to check structure file
        path_in (str): Path to input file
        path_out (str): Path to output file
        
    
    """
    
    command = " ".join(['python', path_to_check_structure, '-i', path_in, '-o', path_out, "--force_save", "--non_interactive", "rem_hydrogen", "--remove Yes"])
    o = os.popen(command).read()
    return o




def select_occupancies(path_to_check_structure, path_in, path_out):
    """
    
    Remove hydrogens from residues
    
    Args:
       
        path_to_check_structure (str): Path to check structure file
        path_in (str): Path to input file
        path_out (str): Path to output file
        
    
    """
    
    command = " ".join(['python', path_to_check_structure, '-i', path_in, '-o', path_out, "--force_save", "--non_interactive", "altloc", "--select occupancy"])
    o = os.popen(command).read()
    return o




def calculate_centroid(arr):
    """
    
    Given an array of 3D points, return the corresponding centroid
    
    Args:
       
        arr (np array/list): Array with 3D points
    
    """
    arr = np.array(arr)
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return np.array([sum_x/length, sum_y/length, sum_z/length])



def extract_ligand_coords(path_to_st, ligand_name, chain):
    """
    
    Extract ligand coordinates.
    
    Args:
       
        path_to_st (str): Path to protein-ligand PDB file
        ligand_name (str): Ligand name (3 letter PDB code)
        chain (str): Ligand chain
        
    CAUTION: If your ligand of interest is repeated within the specified chain, this function sholud not be used.
        
    
    """
    
    structure = PDBParser(QUIET=True).get_structure("st", path_to_st)
    coords = np.array([i.get_coord() for i in structure.get_atoms() if i.get_parent().get_resname() == 'BZU' and i.get_parent().get_parent().id == chain])
    return coords


def create_pocket_centroid(centroid, outfile):
    """
    
    Given a 3D coordinate, create pocket centroid (SD format)
    
    Args:
       
        centroid (3D array): 3D coordinates for pocket centroid (IMPORTANT, SD FORMAT)
        outfile (str): Path to output
        
    CAUTION: Be cautious with PDB file format. Do not change the number of spaces between coordinates.
        
    
    """

    x, y, z = str(round(centroid[0], 3)), str(round(centroid[1], 3)), str(round(centroid[2], 3))
    ctr = " "*(8-len(x)) + x + " "*(8-len(y)) + y + " "*(8-len(z)) + z
    text = """HEADER\nHETATM    1   C  CTR A   1    """ + ctr + """  1.00  1.00           C\nEND"""
    
    # Create PDB file
    with open(outfile.split(".sd")[0] + '.pdb', "w") as f:
        f.write(text)
    
    # Change format (obabel) and remove 
    command = 'obabel ' + outfile.split(".sd")[0] + '.pdb' + " -O " + outfile
    os.system(command)
    os.remove(outfile.split(".sd")[0] + '.pdb')
    
    
    
def create_fpocket_ds(pdb_code, path_to_file):
    """
    
    Create file to rescore fpocket pockets using Prank
    
    Args:
       
        pdb_code (str): PDB code
        path_to_file (str): Path to fpocket ds file
    
    
    """    
    
    
    text = """# Dataset for rescoring Fpocket 3 predictions.
    # (predictions generated with Fpocket 3.1.2 on "./clean" structures)
    #
    # Note: this dataset (unlike fpocket.ds) cannot be used for eval-rescore
    # because protein column points to proteins without known ligands.

    PARAM.PREDICTION_METHOD=fpocket


    HEADER: prediction protein

    """ + pdb_code + """_prepared_out/""" + pdb_code + """_prepared_out.pdb  """ + pdb_code + """_prepared.pdb"""
    with open(path_to_file, "w") as f:
        f.write(text)
        
        
        
        
def read_fpocket_scores(path_to_file):
    
    """
    
    Collect fpocket pocket scores.
    
    Args:
       
        path_to_file (str): Path to fpocket scores  
    
    
    """
    
    file = open(path_to_file, "r").readlines()
    
    fpocket = {}
    for l in file:
        if "Pocket " in l:
            pocket = int(l.split(":")[0].split("Pocket")[1].strip())
            fpocket[pocket] = dict()
        elif l != '\n':
            l = l.split(":")
            feature = l[0].strip()
            value = float(l[1].strip())
            fpocket[pocket][feature] = value
    
    return fpocket




def read_fpocket_centers(path_to_file):
    """
    
    Given a file path (fpocket output pockets pqr), returns it as a python dict
    mapping fpocket pockets with their corresponding alpha sphere centers.
    
    Args:
       
        path_to_file (str): Path to fpocket scores    
    
    """

    centers = {}
    with open(path_to_file, "r") as f:
        for l in f:
            if l.startswith("ATOM"):
                # l = l.split()
                # numb = int(l[4])
                # x = float(l[5])
                # y = float(l[6])
                # z = float(l[7])
                numb, x, y, z = int(l[22:26].strip()), float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip())
                if numb not in centers:
                    centers[numb] = []
                centers[numb].append([x, y, z])
    
    for numb in centers:
        centers[numb] = np.array(centers[numb])
        
    return centers

    
    
def read_prank_scores(path_to_file):
    """
    
    Collect prank pocket scores.
    
    Args:
       
        path_to_file (str): Path to prank scores  
    
    
    """
    
    # Read prank results
    prank = pd.read_csv(path_to_file)
    prank = {i: j for i, j in zip(prank['old_rank'], prank['score'])}
    
    return prank
    
    


    
